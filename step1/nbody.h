/**
 * @file      nbody.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     PCG Assignment 2
 *            N-Body simulation in ACC
 *
 * @version   2021
 *
 * @date      11 November  2020, 11:22 (created) \n
 * @date      11 November  2020, 11:37 (revised) \n
 *
 */

#ifndef __NBODY_H__
#define __NBODY_H__

#include <cstdlib>
#include <cstdio>
#include  <cmath>
#include "h5Helper.h"

/// Gravity constant
constexpr float G = 6.67384e-11f;

/// Collision distance threshold
constexpr float COLLISION_DISTANCE = 0.01f;

/**
 * @struct float4
 * Structure that mimics CUDA float4
 */
struct float4
{
  float x;
  float y;
  float z;
  float w;
};


/// Define sqrtf from CUDA libm library
#pragma acc routine(sqrtf) seq

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Structure with particle data
 */
struct Particles
{
  /// Fill the structure holding the particle/s data
  /// It is recommended to implement constructor / destructor and copyToGPU and copyToCPU routines
  size_t n;       // Number of particles.
  float *pos_x;   // Positions in X
  float *pos_y;   // Positions in Y
  float *pos_z;   // Positions in Z
  float *vel_x;   // Velocities in X
  float *vel_y;   // Velocities in Y
  float *vel_z;   // Velocities in Z
  float *weights; // Weights

  /// Default constructor is not allowed - the particles number must be known in advance
  Particles() = delete;
  /// Copy constructor is not allowed - performance loss due to data copies.
  Particles(const Particles& stc) = delete;
  /// Assignment operator is not allowed - performance loss due to data copies.
  Particles& operator=(const Particles& stc) = delete;

  /**
  * Constructor.
  * This routine construct the particles on both, the CPU and GPU.
  * @param [in] n - Number of particles.
  */
  Particles(size_t n) : n(n)
  {
    /// Construct the particles on the CPU.
    pos_x = new float[n];   // Allocates positions X.
    pos_y = new float[n];   // Allocates positions Y.
    pos_z = new float[n];   // Allocates positions Z.
    vel_x = new float[n];   // Allocates velocities X.
    vel_y = new float[n];   // Allocates velocities Y.
    vel_z = new float[n];   // Allocates velocities Z.
    weights = new float[n]; // Allocates weights.

    /// Construct the particles on the GPU.
    /// Copy the struct into device memory.
    #pragma acc enter data copyin(this)
    /// Allocate the dynamic members into device memory.
    #pragma acc enter data create(pos_x[n])   // Allocates positions X.
    #pragma acc enter data create(pos_y[n])   // Allocates positions Y.
    #pragma acc enter data create(pos_z[n])   // Allocates positions Z.
    #pragma acc enter data create(vel_x[n])   // Allocates velocities X.
    #pragma acc enter data create(vel_y[n])   // Allocates velocities Y.
    #pragma acc enter data create(vel_z[n])   // Allocates velocities Z.
    #pragma acc enter data create(weights[n]) // Allocates weights.
  }/// end of constructor
  ///-----------------------------------------------------------------------------------------------------------

  /**
  * Destructor.
  * This routine destruct the particles on both, the CPU and GPU.
  */
  ~Particles()
  {
    /// First deal with the deallocate of the dynamic members.
    #pragma acc exit data delete(pos_x[n])   // Deallocate positions X.
    #pragma acc exit data delete(pos_y[n])   // Deallocate positions Y.
    #pragma acc exit data delete(pos_z[n])   // Deallocate positions Z.
    #pragma acc exit data delete(vel_x[n])   // Deallocate velocities X.
    #pragma acc exit data delete(vel_y[n])   // Deallocate velocities Y.
    #pragma acc exit data delete(vel_z[n])   // Deallocate velocities Z.
    #pragma acc exit data delete(weights[n]) // Deallocate weights.
    /// Then deallocate the whole struct.
    #pragma acc exit data delete(this)

    /// Destruct the particles on the CPU.
    delete [] pos_x;   // Deallocate positions X.
    delete [] pos_y;   // Deallocate positions Y.
    delete [] pos_z;   // Deallocate positions Z.
    delete [] vel_x;   // Deallocate velocities X.
    delete [] vel_y;   // Deallocate velocities Y.
    delete [] vel_z;   // Deallocate velocities Z.
    delete [] weights; // Deallocate weights.
  } /// Destructor
  ///------------------------------------------------------------------------------------------------------------------

  /// Copy actual CPU particles data to GPU.
  void copyToGPU()
  {
    #pragma acc update device(pos_x[n])   // Copy positions X: CPU => GPU.
    #pragma acc update device(pos_y[n])   // Copy positions Y: CPU => GPU.
    #pragma acc update device(pos_z[n])   // Copy positions Z: CPU => GPU.
    #pragma acc update device(vel_x[n])   // Copy velocities X: CPU => GPU.
    #pragma acc update device(vel_y[n])   // Copy velocities Y: CPU => GPU.
    #pragma acc update device(vel_z[n])   // Copy velocities Z: CPU => GPU.
    #pragma acc update device(weights[n]) // Copy weights: CPU => GPU.
  }/// copyToGPU
  ///------------------------------------------------------------------------------------------------------------------

  /// Copy actual particles data from GPU to CPU.
  void copyToCPU()
  {
    #pragma acc update host(pos_x[n])   // Copy positions X: GPU => CPU.
    #pragma acc update host(pos_y[n])   // Copy positions Y: GPU => CPU.
    #pragma acc update host(pos_z[n])   // Copy positions Z: GPU => CPU.
    #pragma acc update host(vel_x[n])   // Copy velocities X: GPU => CPU.
    #pragma acc update host(vel_y[n])   // Copy velocities Y: GPU => CPU.
    #pragma acc update host(vel_z[n])   // Copy velocities Z: GPU => CPU.
    #pragma acc update host(weights[n]) // Copy weights: GPU => CPU.
  }/// copyToCPU
  ///------------------------------------------------------------------------------------------------------------------

};/// end of Particles
///--------------------------------------------------------------------------------------------------------------------

/**
 * @struct Velocities
 * Velocities of the particles
 */
struct Velocities
{
  /// Fill the structure holding the particle/s data
  /// It is recommended to implement constructor / destructor and copyToGPU and copyToCPU routines
  size_t n;   // Number of particles.
  float *x;   // Velocities X.
  float *y;   // Velocities Y.
  float *z;   // Velocities Z.

  /// Default constructor is not allowed - the velocities number must be known in advance
  Velocities() = delete;
  /// Copy constructor is not allowed - performance loss due to data copies.
  Velocities(const Velocities& stc) = delete;
  /// Assignment operator is not allowed - performance loss due to data copies.
  Velocities& operator=(const Velocities& stc) = delete;

  /**
  * Constructor.
  * This routine construct the velocities on both, the CPU and GPU.
  * @param [in] n - Number of velocities.
  */
  Velocities(size_t n) : n(n)
  {
    /// Construct the velocities on the GPU.
    /// Copy the struct into device memory.
    #pragma acc enter data copyin(this)
    /// Allocate the dynamic members into device memory.
    #pragma acc enter data create(x[n])   // Allocates velocities X.
    #pragma acc enter data create(y[n])   // Allocates velocities Y.
    #pragma acc enter data create(z[n])   // Allocates velocities Z.
  }/// end of constructor
  ///-----------------------------------------------------------------------------------------------------------

  /**
  * Destructor.
  * This routine destruct the velocities on both, the CPU and GPU.
  */
  ~Velocities()
  {
    /// First deal with the deallocate of the dynamic members.
    #pragma acc exit data delete(x[n])   // Deallocate velocities X.
    #pragma acc exit data delete(y[n])   // Deallocate velocities Y.
    #pragma acc exit data delete(z[n])   // Deallocate velocities Z.
    /// Then deallocate the whole struct.
    #pragma acc exit data delete(this)
  } /// Destructor
  ///------------------------------------------------------------------------------------------------------------------

};/// end of Velocities
///----------------------------------------------------------------------------------------------------------------------

/**
 * Compute gravitation velocity
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void calculate_gravitation_velocity(const Particles& p,
                                    Velocities&      tmp_vel,
                                    const int        N,
                                    const float      dt);

/**
 * Calculate collision velocity
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void calculate_collision_velocity(const Particles& p,
                                  Velocities&      tmp_vel,
                                  const int        N,
                                  const float      dt);

/**
 * Update particle position
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void update_particle(const Particles& p,
                     Velocities&      tmp_vel,
                     const int        N,
                     const float      dt);



/**
 * Compute center of gravity - implement in steps 3 and 4.
 * @param [in] p - Particles
 * @param [in] N - Number of particles
 * @return Center of Mass [x, y, z] and total weight[w]
 */
float4 centerOfMassGPU(const Particles& p,
                       const int        N);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Compute center of mass on CPU
 * @param memDesc
 * @return centre of gravity
 */
float4 centerOfMassCPU(MemDesc& memDesc);

#endif /* __NBODY_H__ */
