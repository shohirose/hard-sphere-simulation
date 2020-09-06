#include "shirose/particle_simulator.hpp"

#include <iostream>

#include "shirose/collision.hpp"

namespace shirose {

double EventDrivenParticleSimulator::calcCurrentTime(double oldTime) const {
  const auto currentTimeCandidate = oldTime + deltaTime_;
  const auto currentTime =
      (currentTimeCandidate > endTime_) ? endTime_ : currentTimeCandidate;
  return currentTime;
}

void EventDrivenParticleSimulator::run() {
  using Eigen::Vector2d;
  assert(particles_ && "Null pointer to particles");

  auto& x = particles_->positions;
  auto& v = particles_->velocities;
  const auto& m = particles_->masses;
  const auto numParticles = static_cast<int64_t>(particles_->size());

  bool doUpdateCollisionTime = true;

  // Indices of a colliding particle pair
  std::ptrdiff_t i = 0;
  std::ptrdiff_t j = 0;

  double currentTime = startTime_;
  double oldTime = startTime_;
  double deltaTime = currentTime - oldTime;
  double collisionTime = std::numeric_limits<double>::max();

  // Time loop
  do {
    currentTime = calcCurrentTime(oldTime);
    deltaTime = currentTime - oldTime;

    std::cout << "Current time: " << currentTime << " sec\n";

    // Updates collision time only when required.
    //
    // Because the calculation cost of collision time is expensive, we only
    // recalculate collision time when required. This is when a computed
    // collision time falls in the window of the current time step, i.e.
    // oldTime < collisionTime <= currentTime
    if (doUpdateCollisionTime) {
      std::cout << "Computing time to collision\n";
      const auto timeToCollisionMatrix = calcTimeToCollision(*particles_);
      const auto timeToCollision = timeToCollisionMatrix.minCoeff(&i, &j);
      collisionTime = oldTime + timeToCollision;
    }

    // Checks if the collision time falls in the current time-step window.
    if (collisionTime <= currentTime) {
      std::cout << "Find a collision. Collision time: " << collisionTime
                << " sec\n";

      doUpdateCollisionTime = true;

      // Move current time back to the collision time
      currentTime = collisionTime;
      deltaTime = currentTime - oldTime;

      // Updates positions to the collision time
#ifdef _OPENMP
#pragma omp parallel for
#endif  // _OPENMP
      for (int64_t k = 0; k < numParticles; ++k) {
        x[k] += v[k] * deltaTime;
      }

      // Updates the velocities of colliding particles
      const auto [dv1, dv2] =
          calcVelocityChanges(v[i], v[j], m[i], m[j], restitutionCoeff_);
      v[i] += dv1;
      v[j] += dv2;
    } else {
      doUpdateCollisionTime = false;
    }

    // Updates time
    oldTime = currentTime;

  } while (currentTime < endTime_);

  std::cout << "Simulation has completed." << std::endl;
}

}  // namespace shirose
