#include "shirose/collision.hpp"

#include <limits>

namespace shirose {

namespace {

/// A quadratic equation of the form
/// \f[ a x^2 + 2b x + c = 0 \f]
struct QuadraticEquation {
  double a;
  double b;
  double c;
};

inline double square(double x) noexcept { return x * x; }

inline double calcDeterminant(const QuadraticEquation& eq) noexcept {
  return square(eq.b) - eq.a * eq.c;
}

QuadraticEquation convertToQuadraticEquation(
    const Eigen::Ref<const Eigen::Vector2d>& x1,
    const Eigen::Ref<const Eigen::Vector2d>& x2,
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2,  //
    double r1, double r2) {
  using Eigen::Vector2d;
  const Vector2d dx = x1 - x2;
  const Vector2d dv = v1 - v2;

  const auto a = dv.squaredNorm();
  const auto b = dx.dot(dv);
  const auto c = dx.squaredNorm() - square(r1 + r2);
  return {a, b, c};
}

}  // anonymous namespace

bool doCollide(const Eigen::Ref<const Eigen::Vector2d>& x1,
               const Eigen::Ref<const Eigen::Vector2d>& x2,
               const Eigen::Ref<const Eigen::Vector2d>& v1,
               const Eigen::Ref<const Eigen::Vector2d>& v2,  //
               double r1, double r2) {
  const auto eq = convertToQuadraticEquation(x1, x2, v1, v2, r1, r2);
  const auto det = calcDeterminant(eq);
  return (det >= 0) && (-eq.b - std::sqrt(det) > 0);
}

std::optional<double> calcCollisionTime(
    const Eigen::Ref<const Eigen::Vector2d>& x1,
    const Eigen::Ref<const Eigen::Vector2d>& x2,
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2,  //
    double r1, double r2) {
  const auto eq = convertToQuadraticEquation(x1, x2, v1, v2, r1, r2);
  const auto det = calcDeterminant(eq);

  if (det < 0) {
    return std::nullopt;
  }

  const auto collisionTime = (-eq.b - std::sqrt(det)) / eq.a;

  if (collisionTime > 0) {
    return collisionTime;
  } else {
    return std::nullopt;
  }
}

std::pair<NotAlignedVector2d, NotAlignedVector2d> calcVelocitiesAfterCollision(
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2, double m1, double m2,
    double e) {
  const Eigen::Vector2d dv = v1 - v2;
  const NotAlignedVector2d v1new = v1 - (1 + e) * m2 / (m1 + m2) * dv;
  const NotAlignedVector2d v2new = v2 + (1 + e) * m1 / (m1 + m2) * dv;
  return {v1new, v2new};
}

std::vector<double> calcCollisionTime(const ParticleSystem& particles) {
  const auto numberOfParticles = particles.numberOfParticles();
  std::vector<double> tc(numberOfParticles, std::numeric_limits<double>::max());
  const auto& x = particles.positions;
  const auto& v = particles.velocities;
  const auto& r = particles.radii;

  for (size_t i = 0; i < numberOfParticles; ++i) {
    const auto& x1 = x[i];
    const auto& v1 = v[i];
    const auto r1 = r[i];
    auto& tc1 = tc[i];

    for (size_t j = i + 1; j < numberOfParticles; ++j) {
      const auto& x2 = x[j];
      const auto& v2 = v[j];
      const auto r2 = r[j];
      auto& tc2 = tc[j];

      const auto tc12opt = calcCollisionTime(x1, x2, v1, v2, r1, r2);
      if (tc12opt) {
        const auto tc12 = tc12opt.value();
        if (tc12 < tc1) tc1 = tc12;
        if (tc12 < tc2) tc2 = tc12;
      }
    }
  }

  return tc;
}

}  // namespace shirose
