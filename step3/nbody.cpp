/**
 * @file      nbody.cpp
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

#include <math.h>
#include <cfloat>
#include "nbody.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute gravitation and collision velocity, and subsequently update particle positions
void calculate_velocity(const Particles &p_in, Particles &p_out, const int N, const float dt) {

  /// Initialise of auxiliary accumulators of velocity
  float tmp_vel_x = 0.0f;
  float tmp_vel_y = 0.0f;
  float tmp_vel_z = 0.0f;

  #pragma acc parallel loop present(p_in, p_out) gang worker vector
  for (unsigned int i = 0; i < N; i++) {
  #pragma acc loop seq
    /// The iterations over all particles to compute the gravitation velocity to them
    for (int j = 0; j < N; j++) {

      /// Instruction Level Parallelism
      float s = -G * dt * p_in.weights[j];
      /// Computes the distance between the relevant particles
      float r_x = p_in.pos_x[i] - p_in.pos_x[j];
      float r_y = p_in.pos_y[i] - p_in.pos_y[j];
      float r_z = p_in.pos_z[i] - p_in.pos_z[j];

      /// Computes inverse distance between particles and their distances
      float inv_dist = sqrtf(r_x * r_x + r_y * r_y + r_z * r_z);
      /// Computes the gravitation velocity (Fg_dt_m2_r)
      s /= (inv_dist * inv_dist * inv_dist + FLT_MIN);

      /// The speed that a particle body receives due to the strength of the relevant particle
      tmp_vel_x += (inv_dist > COLLISION_DISTANCE) ? r_x * s : 0.0f;
      tmp_vel_y += (inv_dist > COLLISION_DISTANCE) ? r_y * s : 0.0f;
      tmp_vel_z += (inv_dist > COLLISION_DISTANCE) ? r_z * s : 0.0f;

      /// Checks whether the particles are in the sufficient near distance for collision
      if (inv_dist > 0.0f && inv_dist < COLLISION_DISTANCE) {
        /// Computes the temporary partial results to eliminate recalculation
        float weight_diff = p_in.weights[i] - p_in.weights[j];
        float weight_sum = p_in.weights[i] + p_in.weights[j];
        float weight_j_x_2 = 2 * p_in.weights[j];

        /// Computes the collision velocities between the relevant particles and accumulate the results
        tmp_vel_x += inv_dist > 0.0f ?
                     ((weight_diff * p_in.vel_x[i] + weight_j_x_2 * p_in.vel_x[j]) / weight_sum) - p_in.vel_x[i] : 0.0f;
        tmp_vel_y += inv_dist > 0.0f ?
                     ((weight_diff * p_in.vel_y[i] + weight_j_x_2 * p_in.vel_y[j]) / weight_sum) - p_in.vel_y[i] : 0.0f;
        tmp_vel_z += inv_dist > 0.0f ?
                     ((weight_diff * p_in.vel_z[i] + weight_j_x_2 * p_in.vel_z[j]) / weight_sum) - p_in.vel_z[i] : 0.0f;
      }
    }

    /// Updates the velocity of particles with respect to the computed gravitation and collision velocity
    p_out.vel_x[i] = p_in.vel_x[i] + tmp_vel_x;
    p_out.vel_y[i] = p_in.vel_y[i] + tmp_vel_y;
    p_out.vel_z[i] = p_in.vel_z[i] + tmp_vel_z;

    /// Updates the positions of particles with respect to the updated velocity
    p_out.pos_x[i] = p_in.pos_x[i] + (p_in.vel_x[i] + tmp_vel_x) * dt;
    p_out.pos_y[i] = p_in.pos_y[i] + (p_in.vel_y[i] + tmp_vel_y) * dt;
    p_out.pos_z[i] = p_in.pos_z[i] + (p_in.vel_z[i] + tmp_vel_z) * dt;

    /// Set the null before the next computation iteration  of the velocity accumulation
    tmp_vel_x = 0.0f;
    tmp_vel_y = 0.0f;
    tmp_vel_z = 0.0f;
  }

}


/// Compute center of gravity
float4 centerOfMassGPU(const Particles &p, const int N) {
  float sum_x = 0.0;  // Initialize overall sum of Positions X.
  float sum_y = 0.0;  // Initialize overall sum of Positions Y.
  float sum_z = 0.0;  // Initialize overall sum of Positions Z.
  float sum_w = 0.0;  // Initialize overall sum of Weights.

  /// Reduction across Multi-level Thread Parallelism in the Same Loop
  #pragma acc parallel loop present(p) reduction(+:sum_x, sum_y, sum_z, sum_w) gang worker vector
  for (unsigned int i = 0; i < N; i++) {
    sum_x += (p.pos_x[i] * p.weights[i]);   // Add the weighting position X of the current particle.
    sum_y += (p.pos_y[i] * p.weights[i]);   // Add the weighting position Y of the current particle.
    sum_z += (p.pos_z[i] * p.weights[i]);   // Add the weighting position Z of the current particle.
    sum_w += p.weights[i];                  // Add the weights of the current particle.
  }

  /// Compute the final results of Center of Mass
  /// x_com = (x1*w1 + x2*w2 + ... xn+wn) / (w1 + w2 + ... + wn)
  return {sum_x / sum_w, sum_y / sum_w, sum_z / sum_w, sum_w};
}/// end of centerOfMassGPU
///---------------------------------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute center of mass on CPU
float4 centerOfMassCPU(MemDesc &memDesc) {
  float4 com = {0, 0, 0, 0};

  for (int i = 0; i < memDesc.getDataSize(); i++) {
    // Calculate the vector on the line connecting points and most recent position of center-of-mass
    const float dx = memDesc.getPosX(i) - com.x;
    const float dy = memDesc.getPosY(i) - com.y;
    const float dz = memDesc.getPosZ(i) - com.z;

    // Calculate weight ratio only if at least one particle isn't massless
    const float dw = ((memDesc.getWeight(i) + com.w) > 0.0f)
                     ? (memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w)) : 0.0f;

    // Update position and weight of the center-of-mass according to the weight ration and vector
    com.x += dx * dw;
    com.y += dy * dw;
    com.z += dz * dw;
    com.w += memDesc.getWeight(i);
  }
  return com;
}// end of centerOfMassCPU
//----------------------------------------------------------------------------------------------------------------------
