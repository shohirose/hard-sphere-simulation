#ifndef SHIROSE_PARTICLE_SIMULATOR_HPP
#define SHIROSE_PARTICLE_SIMULATOR_HPP

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <tuple>

#include "shirose/particle_system.hpp"

namespace shirose {

class EventDrivenParticleSimulator {
 public:
  EventDrivenParticleSimulator() = default;

  void setTime(double startTime, double endTime, double deltaTime) noexcept {
    startTime_ = startTime;
    endTime_ = endTime;
    deltaTime_ = deltaTime;
  }

  void setRestitutionCoeff(double restitutionCoeff) noexcept {
    restitutionCoeff_ = restitutionCoeff;
  }

  void setParticles(const std::shared_ptr<ParticleSystem>& particles) {
    particles_ = particles;
  }

  void run();

 private:
  double calcCurrentTime(double oldTime) const;

  std::shared_ptr<ParticleSystem> particles_;
  double startTime_;
  double endTime_;
  double deltaTime_;
  double restitutionCoeff_;
};

}  // namespace shirose

#endif  // SHIROSE_PARTICLE_SIMULATOR_HPP