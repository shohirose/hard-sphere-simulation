#ifndef SHIROSE_VECTOR_HPP
#define SHIROSE_VECTOR_HPP

#include <Eigen/Core>

namespace shirose {

using NotAlignedVector2d = Eigen::Matrix<double, 2, 1, Eigen::DontAlign>;

}  // namespace shirose

#endif  // SHIROSE_VECTOR_HPP