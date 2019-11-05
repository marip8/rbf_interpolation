#ifndef RBF_INTERPOLATION_RBF_MATRIX_TYPES_H
#define RBF_INTERPOLATION_RBF_MATRIX_TYPES_H

#include <Eigen/Geometry>

namespace rbf_interpolation
{
template<typename T>
using RBFVectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using RBFMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using DataMap = std::vector < std::pair<const T, RBFVectorX<T>>>;

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_MATRIX_TYPES_H
