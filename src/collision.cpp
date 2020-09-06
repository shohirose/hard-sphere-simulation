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

/// @brief Converts the collision condition of two particles into a quadratic
/// equation.
/// @param[in] x1 Position of particle 1
/// @param[in] x2 Position of particle 2
/// @param[in] v1 Velocity of particle 1
/// @param[in] v2 Velocity of particle 2
/// @param[in] r1 Radius of particle 1
/// @param[in] r2 Radius of particle 2
inline QuadraticEquation convertToQuadraticEquation(
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

constexpr double eps = 1e-10;

}  // anonymous namespace

bool doCollide(const Eigen::Ref<const Eigen::Vector2d>& x1,
               const Eigen::Ref<const Eigen::Vector2d>& x2,
               const Eigen::Ref<const Eigen::Vector2d>& v1,
               const Eigen::Ref<const Eigen::Vector2d>& v2,  //
               double r1, double r2) {
  const auto eq = convertToQuadraticEquation(x1, x2, v1, v2, r1, r2);
  const auto det = calcDeterminant(eq);
  const auto timeToCollision = (-eq.b - std::sqrt(det)) / eq.a;
  return det >= 0 && timeToCollision > eps;
}

std::pair<bool, double> calcTimeToCollision(
    const Eigen::Ref<const Eigen::Vector2d>& x1,
    const Eigen::Ref<const Eigen::Vector2d>& x2,
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2,  //
    double r1, double r2) {
  const auto eq = convertToQuadraticEquation(x1, x2, v1, v2, r1, r2);
  const auto det = calcDeterminant(eq);

  if (det < 0) {
    return {false, 0};
  }

  const auto timeToCollision = (-eq.b - std::sqrt(det)) / eq.a;

  if (timeToCollision > eps) {
    return {true, timeToCollision};
  } else {
    return {false, 0};
  }
}

Pair<Eigen::Vector2d, Eigen::Vector2d> calcVelocityChanges(
    const Eigen::Ref<const Eigen::Vector2d>& v1,
    const Eigen::Ref<const Eigen::Vector2d>& v2, double m1, double m2,
    double e) {
  using Eigen::Vector2d;
  const Vector2d dv = v1 - v2;
  const Vector2d dv1 = -(1 + e) * m2 / (m1 + m2) * dv;
  const Vector2d dv2 = (1 + e) * m1 / (m1 + m2) * dv;
  return {dv1, dv2};
}

Eigen::MatrixXd calcTimeToCollision(const ParticleSystem& particles) {
  const auto& x = particles.positions;
  const auto& v = particles.velocities;
  const auto& r = particles.radii;

  using Eigen::MatrixXd;
  const auto rows = static_cast<int64_t>(particles.size());
  const auto cols = rows;
  // Collision time matrix
  MatrixXd tc =
      MatrixXd::Constant(rows, cols, std::numeric_limits<double>::max());

#ifdef _OPENMP
#pragma omp parallel for
#endif  // _OPENMP
  for (int64_t j = 0; j < cols; ++j) {
    const auto& x1 = x[j];
    const auto& v1 = v[j];
    const auto r1 = r[j];

    for (int64_t i = j + 1; i < rows; ++i) {
      const auto& x2 = x[i];
      const auto& v2 = v[i];
      const auto r2 = r[i];

      const auto [doCollide, timeToCollision] =
          calcTimeToCollision(x1, x2, v1, v2, r1, r2);

      if (doCollide && (timeToCollision < tc(i, j))) {
        tc(i, j) = timeToCollision;
        tc(j, i) = timeToCollision;
      }
    }
  }

  return tc;
}

}  // namespace shirose
