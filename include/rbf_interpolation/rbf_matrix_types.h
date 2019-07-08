#ifndef RBF_INTERPOLATION_RBF_MATRIX_TYPES_H
#define RBF_INTERPOLATION_RBF_MATRIX_TYPES_H

#include <Eigen/Eigen>

namespace rbf_interpolation
{

template<typename T>
struct RBFVectorX
{
  typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> type;
};

template<typename T>
struct RBFMatrixX
{
  typedef typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> type;
};

template<typename T>
struct DataMap
{
  typedef typename std::vector<std::pair<const T, typename RBFVectorX<T>::type>> type;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_MATRIX_TYPES_H
