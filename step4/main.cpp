/**
 * @file      main.cpp
 *
 * @author    Simon Stupinsky; xstupi00 \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xstupi00@fit.vutbr.cz
 *
 * @brief     PCG Assignment 2
 *            N-Body simulation in ACC
 *
 * @version   2021
 *
 * @date      11 November  2020, 11:22 (created) \n
 * @date      21 December  2020, 21:10 (revised) \n
 *
 */

#include <chrono>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "nbody.h"
#include "h5Helper.h"

/**
 * Main routine of the project
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  // Parse command line parameters
  if (argc != 7) {
    printf("Usage: nbody <N> <dt> <steps> <write intesity> <input> <output>\n");
    exit(EXIT_FAILURE);
  }

  const int N = std::stoi(argv[1]);
  const float dt = std::stof(argv[2]);
  const int steps = std::stoi(argv[3]);
  const int writeFreq = (std::stoi(argv[4]) > 0) ? std::stoi(argv[4]) : 0;

  printf("N: %d\n", N);
  printf("dt: %f\n", dt);
  printf("steps: %d\n", steps);

  const size_t recordsNum = (writeFreq > 0) ? (steps + writeFreq - 1) / writeFreq : 0;


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                         Code to be implemented                                                   //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // 1.  Memory allocation on CPU
  Particles particles_in(N);
  Particles particles_out(N);

  // 2. Create memory descriptor
  /*
   * Caution! Create only after CPU side allocation
   * parameters:
   *                                    Stride of two               Offset of the first
   *             Data pointer           consecutive elements        element in floats,
   *                                    in floats, not bytes        not bytes
  */
  MemDesc md(
    particles_in.pos_x, 1, 0,            // Position in X
    particles_in.pos_y, 1, 0,            // Position in Y
    particles_in.pos_z, 1, 0,            // Position in Z
    particles_in.vel_x, 1, 0,            // Velocity in X
    particles_in.vel_y, 1, 0,            // Velocity in Y
    particles_in.vel_z, 1, 0,            // Velocity in Z
    particles_in.weights, 1, 0,            // Weights
    N,                                                                           // Number of particles
    recordsNum                                                                   // Number of records in output file
  );

  H5Helper h5Helper(argv[5], argv[6], md);

  // Read data
  try {
    h5Helper.init();
    h5Helper.readParticleData();
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // 3. Copy data to GPU
  particles_in.copyToGPU();

  // Copy weights from input to output particles on CPU due to the first iteration of the computation.
  memcpy(particles_out.weights, particles_in.weights, N * sizeof(float));
  // Copy weights to output particles on GPU.
  #pragma acc update device(particles_out.weights[N]) // Copy weights: CPU => GPU.

  // Alloc the structure for Center of Mass data on the CPU.
  float4 *__restrict comOnGPU = (float4 *) malloc(sizeof(float4));
  // Alloc the structure for Center of Mass data on the GPU.
  #pragma acc enter data create(comOnGPU[0])

  // Start the time
  auto startTime = std::chrono::high_resolution_clock::now();

  // 4. Run the loop - calculate new Particle positions.
  for (int s = 0; s < steps; s++) {

    /// Computes Center of Mass from the particles data computed in the previous iteration.
    centerOfMassGPU((s & 1) ? particles_out : particles_in, comOnGPU, N);

    /// In step 4 - fill in the code to store Particle snapshots.
    if (writeFreq > 0 && (s % writeFreq == 0)) {

      if (s > 0) {
        /// Copy the particles data from GPU to CPU computed in the previous iteration.
        (s & 1ul) ? particles_out.copyToCPUAsync() : particles_in.copyToCPUAsync();
      }

      /// Wait for the finish of COM data computation within the current iteration.
      #pragma acc update host(comOnGPU[0]) async(UPDATE_COM_STREAM) wait(SERIAL_COM_STREAM)

      /// Wait for the particle data update from GPU to CPU.
      #pragma acc wait(UPDATE_PARTICLES_STREAM)

      /// h5helper reads data only from the input particle data structure
      /// Therefore, we need to copy output particle when they was as input in the previous iteration
      if (s & 1ul) {
        std::memcpy(particles_in.pos_x, particles_out.pos_x, sizeof(float) * N);
        std::memcpy(particles_in.pos_y, particles_out.pos_y, sizeof(float) * N);
        std::memcpy(particles_in.pos_z, particles_out.pos_z, sizeof(float) * N);
        std::memcpy(particles_in.vel_x, particles_out.vel_x, sizeof(float) * N);
        std::memcpy(particles_in.vel_y, particles_out.vel_y, sizeof(float) * N);
        std::memcpy(particles_in.vel_z, particles_out.vel_z, sizeof(float) * N);
        std::memcpy(particles_in.weights, particles_out.weights, sizeof(float) * N);
      }

      /// Writing particles data to the file
      h5Helper.writeParticleData(s / writeFreq);

      /// Wait for the COM data update from GPU to CPU
      #pragma acc wait(UPDATE_COM_STREAM)
      /// Writing computed values Center of Mass to the file
      h5Helper.writeCom(comOnGPU[0].x, comOnGPU[0].y, comOnGPU[0].z, comOnGPU[0].w, s / writeFreq);
    }

    /// Compute gravitation and collision velocity, and subsequently update particle positions
    calculate_velocity(
      (s & 1ul) ? particles_out : particles_in,
      (s & 1ul) ? particles_in : particles_out,
      N, dt
    );

    /// Synchronization of the both compute streams across the iterations of the loop.
    /// Both compute streams have to wait for each other finish the current iteration.
    #pragma acc wait(COMPUTE_VELOCITIES_STREAM) async(COMPUTE_COM_STREAM)
    #pragma acc wait(COMPUTE_COM_STREAM) async(COMPUTE_VELOCITIES_STREAM)
  }// for s ...
  // Wait for finish of all computations within async streams
  #pragma acc wait

  // 5. In steps 3 and 4 -  Compute center of gravity
  centerOfMassGPU((steps & 1) ? particles_out : particles_in, comOnGPU, N);

  // Stop watchclock
  const auto endTime = std::chrono::high_resolution_clock::now();
  const double time = (endTime - startTime) / std::chrono::milliseconds(1);
  printf("Time: %f s\n", time / 1000);

  // 5. Copy data from GPU back to CPU.
  (steps & 1) ? particles_out.copyToCPU() : particles_in.copyToCPU(); // Copy particles data from GPU to CPU.
  #pragma acc update host(comOnGPU[0]) wait(SERIAL_COM_STREAM)        // Copy COM data from GPU to CPU.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Calculate center of gravity
  float4 comOnCPU = centerOfMassCPU(md);

  std::cout << "Center of mass on CPU:" << std::endl
            << comOnCPU.x << ", "
            << comOnCPU.y << ", "
            << comOnCPU.z << ", "
            << comOnCPU.w
            << std::endl;

  std::cout << "Center of mass on GPU:" << std::endl
            << comOnGPU[0].x << ", "
            << comOnGPU[0].y << ", "
            << comOnGPU[0].z << ", "
            << comOnGPU[0].w
            << std::endl;

  // Store final positions of the particles into a file
  h5Helper.writeComFinal(comOnGPU[0].x, comOnGPU[0].y, comOnGPU[0].z, comOnGPU[0].w);
  h5Helper.writeParticleDataFinal();

  // Deallocate the Center of Mass data structure on the GPU.
  #pragma acc exit data delete(comOnGPU[0])
  // Deallocate the Center of Mass data structure on the CPU.
  std::free(comOnGPU);

  return EXIT_SUCCESS;
}// end of main
//----------------------------------------------------------------------------------------------------------------------

