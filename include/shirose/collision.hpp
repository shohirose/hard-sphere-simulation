#ifndef SHIROSE_COLLISION_HPP
#define SHIROSE_COLLISION_HPP

#include <Eigen/Core>
#include <utility>
#include <vector>

#include "shirose/particle_system.hpp"
#include "shirose/vector.hpp"

namespace shirose {

/// @brief If two particles will collide.
/// @param[in] x1 Position of particle 1
/// @param[in] x2 Position of particle 2
/// @param[in] v1 Velocity of particle 1
/// @param[in] v2 Velocity of particle 2
/// @param[in] r1 Radius of particle 1
/// @param[in] r2 Radius of particle 2
bool doCollide(const Eigen::Ref<const Eigen::Vector2d>& x1,
               const Eigen::Ref<const Eigen::Vector2d>& x2,
               const Eigen::Ref<const Eigen::Vector2d>& v1,
               const Eigen::Ref<const Eigen::Vector2d>& v2,  //
               double r1, double r2);

/// @brief Computes collision time of two particles.
/// @param[in] x1 Position of particle 1
/// @param[in] x2 Position of particle 2
/// @param[in] v1 Velocity of particle 1
/// @param[in] v2 Velocity of particle 2
/// @param[in] r1 Radius of particle 1
/// @param[in] r2 Radius of particle 2
///
/// If two particles will not collide, the maximum value of double is returned.
double calcCollisionTime(const Eigen::Ref<const Eigen::Vector2d>& x1,
                         const Eigen::Ref<const Eigen::Vector2d>& x2,
                         const Eigen::Ref<const Eigen::Vector2d>& v1,
                         const Eigen::Ref<const Eigen::Vector2d>& v2,  //
                         double r1, double r2);

/// @brief Computes collision times between particles.
/// @param[in] particles Particles
///
/// If a particle will not collide, its collision time is set to the maximum
/// value of double.
std::vector<double> calcCollisionTime(const ParticleSystem& particles);

/// @brief Comutes velocities after collision.
/// @param[in] v1 Velocity of particle 1
/// @param[in] v2 Velocity of particle 2
/// @param[in] m1 Mass of particle 1
/// @param[in] m2 Mass of particle 2
/// @param[in] e Coefficient of restitution
std::pair<NotAlignedVector2d, NotAlignedVector2d> calcVelocitiesAfterCollision(
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2,  //
    double m1, double m2, double e);

}  // namespace shirose

#endif  // SHIROSE_COLLISION_HPP
