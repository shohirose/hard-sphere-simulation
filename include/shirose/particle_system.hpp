#ifndef SHIROSE_PARTICLE_SYSTEM_HPP
#define SHIROSE_PARTICLE_SYSTEM_HPP

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <type_traits>
#include <cassert>

namespace shirose {

struct ParticleSystem {
  template <typename T>
  using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;

  /// @brief Positions of particles at the current time step
  StdVector<Eigen::Vector2d> positions;

  /// @brief Velocities of particles at the current time step
  StdVector<Eigen::Vector2d> velocities;

  /// @brief Positions of particles at the previous time step
  StdVector<Eigen::Vector2d> previousPositions;

  /// @brief Velocities of particles at the previous time step
  StdVector<Eigen::Vector2d> previousVelocities;

  /// @brief Masses of particles
  std::vector<double> masses;

  /// @brief Radii of particles
  std::vector<double> radii;

   /// @brief Returns the number of particles.
  std::size_t numberOfParticles() const noexcept {
    assert(this->positions.size() == this->velocities.size() &&
           this->positions.size() == this->masses.size() &&
           this->positions.size() == this->radii.size());
    return this->positions.size();
  }

  /// @brief Stores the current velocities and positions in the previous ones.
  void storePreviousTimeStep() noexcept {
    this->previousPositions = this->positions;
    this->previousVelocities = this->velocities;
  }
};

}  // namespace shirose

#endif  // SHIROSE_PARTICLE_SYSTEM_HPP