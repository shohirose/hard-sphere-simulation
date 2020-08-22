#ifndef SHIROSE_PARTICLE_SYSTEM_HPP
#define SHIROSE_PARTICLE_SYSTEM_HPP

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <type_traits>
#include <cassert>

namespace shirose {

struct ParticleSystem {
  /// Use Eigen::aligned allocator for Eigen::Vector2d.
  /// This is because Eigen::Vector2d requires aligned memory allocator.
  /// Please refer to:
  /// https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
  template <typename T>
  using StdVector = std::vector<T, Eigen::aligned_allocator<T>>;

  StdVector<Eigen::Vector2d> positions;
  StdVector<Eigen::Vector2d> velocities;
  std::vector<double> masses;
  std::vector<double> radii;

   /// @brief Returns the number of particles.
  std::size_t size() const noexcept {
    assert(this->positions.size() == this->velocities.size() &&
           this->positions.size() == this->masses.size() &&
           this->positions.size() == this->radii.size());
    return this->positions.size();
  }
};

}  // namespace shirose

#endif  // SHIROSE_PARTICLE_SYSTEM_HPP