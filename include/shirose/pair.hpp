#ifndef SHIROSE_PAIR_HPP
#define SHIROSE_PAIR_HPP

#include <Eigen/Core>

namespace shirose {

template <typename T1, typename T2>
struct Pair {
  T1 first;
  T2 second;
};

/// @brief Partial specialization for Eigen::Matrix.
///
/// Please refer to:
/// https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
template <typename T1, int Rows1, int Cols1, int Options1, typename T2, int Rows2, int Cols2, int Options2>
struct Pair<Eigen::Matrix<T1, Rows1, Cols1, Options1>, Eigen::Matrix<T2, Rows2, Cols2, Options2>> {
  using Matrix1 = Eigen::Matrix<T1, Rows1, Cols1, Options1>;
  using Matrix2 = Eigen::Matrix<T2, Rows2, Cols2, Options2>;

  static constexpr bool needsToAlign =
      (sizeof(Matrix1) % 16) == 0 || (sizeof(Matrix2) % 16) == 0;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(needsToAlign)

  Matrix1 first;
  Matrix2 second;
};

template <typename T1, typename T2>
Pair<T1, T2> makePair(T1&& first, T2&& second) {
  return {first, second};
}

}  // namespace shirose

#endif  // SHIROSE_PAIR_HPP