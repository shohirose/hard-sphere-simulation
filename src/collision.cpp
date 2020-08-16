#include "shirose/collision.hpp"

#include <limits>

namespace shirose {

bool doCollide(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2,
               const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, double r1,
               double r2) {
  using Eigen::Vector2d;
  Vector2d dx = x1 - x2;
  Vector2d dv = v1 - v2;

  const auto xv = dx.dot(dv);
  const auto vsqr = dv.squaredNorm();
  const auto xsqr = dx.squaredNorm();
  const auto rsum = r1 + r2;
  const auto det = xv * xv - vsqr * (xsqr - rsum * rsum);

  return (det > 0) && (-(xv + std::sqrt(det)) > 0);
}

double calcCollisionTime(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2,
                         const Eigen::Vector2d& v1, const Eigen::Vector2d& v2,
                         double r1, double r2) {
  using Eigen::Vector2d;
  Vector2d dx = x1 - x2;
  Vector2d dv = v1 - v2;

  const auto xv = dx.dot(dv);
  const auto vsqr = dv.squaredNorm();
  const auto xsqr = dx.squaredNorm();
  const auto rsum = r1 + r2;
  const auto det = xv * xv - vsqr * (xsqr - rsum * rsum);

  return (det > 0) ? -(xv + std::sqrt(det)) / vsqr
                   : std::numeric_limits<double>::max();
}

std::pair<NotAlignedVector2d, NotAlignedVector2d> calcVelocitiesAfterCollision(
    const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, double m1, double m2,
    double e) {
  const Eigen::Vector2d dv = v1 - v2;
  const NotAlignedVector2d v1new = v1 - (1 + e) * m2 / (m1 + m2) * dv;
  const NotAlignedVector2d v2new = v2 + (1 + e) * m1 / (m1 + m2) * dv;
  return {v1new, v2new};
}

std::vector<double> calcCollisionTime(const ParticleSystem& particles) {
  const auto numberOfParticles =
      static_cast<int64_t>(particles.numberOfParticles());
  std::vector<double> tc(numberOfParticles, std::numeric_limits<double>::max());
  const auto& x = particles.positions;
  const auto& v = particles.velocities;
  const auto& r = particles.radii;

  for (int64_t i = 0; i < numberOfParticles; ++i) {
    const auto& x1 = x[i];
    const auto& v1 = v[i];
    const auto r1 = r[i];
    auto& tc1 = tc[i];

    for (int64_t j = i + 1; j < numberOfParticles; ++j) {
      const auto& x2 = x[j];
      const auto& v2 = v[j];
      const auto r2 = r[j];
      auto& tc2 = tc[j];

      const auto tc12 = calcCollisionTime(x1, x2, v1, v2, r1, r2);
      if (tc12 < tc1) tc1 = tc12;
      if (tc12 < tc2) tc2 = tc12;
    }
  }
}

}  // namespace shirose
