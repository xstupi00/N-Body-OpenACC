/**
 * @file      nbody.cpp
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
 * @date      21 December  2020, 17:56 (revised) \n
 *
 */

#include <math.h>
#include <cfloat>
#include "nbody.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute gravitation velocity
void calculate_gravitation_velocity(
  const Particles &p, Velocities &tmp_vel, const int N, const float dt
) {

  #pragma acc parallel loop present(p, tmp_vel) gang worker vector
  for (unsigned int i = 0; i < N; i++) {

    /// Initialise of auxiliary accumulators of velocity
    float tmp_vel_x = 0.0f;
    float tmp_vel_y = 0.0f;
    float tmp_vel_z = 0.0f;

  #pragma acc loop seq
    for (unsigned int j = 0; j < N; j++) {

      /// ILP
      float s = -G * dt * p.weights[j];
      /// Computes the distance between the relevant particles
      float r_x = p.pos_x[i] - p.pos_x[j];
      float r_y = p.pos_y[i] - p.pos_y[j];
      float r_z = p.pos_z[i] - p.pos_z[j];

      /// Computes inverse distance between particles and their distances
      float inv_dist = sqrtf(r_x * r_x + r_y * r_y + r_z * r_z);
      /// Computes the gravitation velocity (Fg_dt_m2_r)
      s /= (inv_dist * inv_dist * inv_dist + FLT_MIN);

      /// The speed that a particle body receives due to the strength of the relevant particle
      tmp_vel_x += (inv_dist > COLLISION_DISTANCE) ? r_x * s : 0.0f;
      tmp_vel_y += (inv_dist > COLLISION_DISTANCE) ? r_y * s : 0.0f;
      tmp_vel_z += (inv_dist > COLLISION_DISTANCE) ? r_z * s : 0.0f;
    }

    /// Stores the final computed velocities of each particle to the auxiliary velocity vector
    tmp_vel.x[i] = tmp_vel_x;
    tmp_vel.y[i] = tmp_vel_y;
    tmp_vel.z[i] = tmp_vel_z;
  }

}/// end of calculate_gravitation_velocity
///---------------------------------------------------------------------------------------------------------------------

/// Compute collision velocity
void calculate_collision_velocity(
  const Particles &p, Velocities &tmp_vel, const int N, const float dt
) {

  #pragma acc parallel loop present(p, tmp_vel) gang worker vector
  for (unsigned int i = 0; i < N; i++) {

    /// Initialise of auxiliary accumulators of velocity
    float tmp_vel_x = 0.0f;
    float tmp_vel_y = 0.0f;
    float tmp_vel_z = 0.0f;

  #pragma acc loop seq
    for (unsigned int j = 0; j < N; j++) {

      /// Computes the distance between the relevant particles
      float r_x = p.pos_x[i] - p.pos_x[j];
      float r_y = p.pos_y[i] - p.pos_y[j];
      float r_z = p.pos_z[i] - p.pos_z[j];

      /// Computes inverse distance between particles and their distances
      float inv_dist = sqrtf(r_x * r_x + r_y * r_y + r_z * r_z);

      /// Checks whether the particles are in the sufficient near distance for collision
      if (inv_dist > 0.0f && inv_dist < COLLISION_DISTANCE) {
        /// Computes the temporary partial results to eliminate recalculation
        float weight_diff = p.weights[i] - p.weights[j];
        float weight_sum = p.weights[i] + p.weights[j];
        float weight_j_x_2 = 2 * p.weights[j];

        /// Computes the collision velocities between the relevant particles and accumulate the results
        tmp_vel_x += inv_dist > 0.0f ?
                     ((weight_diff * p.vel_x[i] + weight_j_x_2 * p.vel_x[j]) / weight_sum) - p.vel_x[i] : 0.0f;
        tmp_vel_y += inv_dist > 0.0f ?
                     ((weight_diff * p.vel_y[i] + weight_j_x_2 * p.vel_y[j]) / weight_sum) - p.vel_y[i] : 0.0f;
        tmp_vel_z += inv_dist > 0.0f ?
                     ((weight_diff * p.vel_z[i] + weight_j_x_2 * p.vel_z[j]) / weight_sum) - p.vel_z[i] : 0.0f;
      }
    }

    /// Stores the final computed velocities of each particle to the auxiliary velocity vector
    tmp_vel.x[i] += tmp_vel_x;
    tmp_vel.y[i] += tmp_vel_y;
    tmp_vel.z[i] += tmp_vel_z;
  }

}/// end of calculate_collision_velocity
///---------------------------------------------------------------------------------------------------------------------

/// Update particle position
void update_particle(
  const Particles &p, Velocities &tmp_vel, const int N, const float dt
) {

  #pragma acc parallel loop present(p, tmp_vel)
  for (unsigned int i = 0; i < N; i++) {
    /// Updates the velocity of particles with respect to the computed gravitation and collision velocity
    p.vel_x[i] += tmp_vel.x[i];
    p.vel_y[i] += tmp_vel.y[i];
    p.vel_z[i] += tmp_vel.z[i];

    /// Updates the positions of particles with respect to the updated velocity
    p.pos_x[i] += p.vel_x[i] * dt;
    p.pos_y[i] += p.vel_y[i] * dt;
    p.pos_z[i] += p.vel_z[i] * dt;
  }

}/// end of update_particle
///---------------------------------------------------------------------------------------------------------------------



/// Compute center of gravity
float4 centerOfMassGPU(const Particles &p,
                       const int N) {

  return {0.0f, 0.0f, 0.0f, 0.0f};
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
