#include "shirose/collision.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

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
  const auto& x = particles.positions;
  const auto& v = particles.velocities;
  const auto& r = particles.radii;

  using Eigen::MatrixXd;
  const auto numberOfParticles = static_cast<int64_t>(particles.numberOfParticles());
  // Collision time matrix
  MatrixXd tcMatrix = MatrixXd::Constant(numberOfParticles, numberOfParticles,
                                         std::numeric_limits<double>::max());

  #ifdef _OPENMP
#pragma omp parallel for
  #endif // _OPENMP
  for (int64_t j = 0; j < numberOfParticles; ++j) {
    const auto& x1 = x[j];
    const auto& v1 = v[j];
    const auto r1 = r[j];

    for (int64_t i = j + 1; i < numberOfParticles; ++i) {
      const auto& x2 = x[i];
      const auto& v2 = v[i];
      const auto r2 = r[i];

      const auto doCollide = calcCollisionTime(x1, x2, v1, v2, r1, r2);
      if (doCollide) {
        const auto tc12 = doCollide.value();
        if (tc12 < tcMatrix(i, j)) {
          tcMatrix(i, j) = tc12;
          tcMatrix(j, i) = tc12;
        }
      }
    }
  }

  // Collision time of each particle
  std::vector<double> tc(numberOfParticles);

  Eigen::Map<Eigen::VectorXd> tcMap(tc.data(), numberOfParticles);

  // Collision time of each particle is the minimum value of each column of the
  // collision time matrix.
  tcMap = tcMatrix.colwise().minCoeff().transpose();

  return tc;
}

}  // namespace shirose
