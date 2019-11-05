#ifndef RBF_INTERPOLATION_RBF_SOLVER_H
#define RBF_INTERPOLATION_RBF_SOLVER_H

#include <map>
#include <functional>
#include <memory>
#include "rbf_interpolation/rbf_base.h"

namespace rbf_interpolation
{

template <typename T>
class RBFSolver
{
public:
  RBFSolver(const typename RBFBase<T>::Ptr rbf);

  /**
   * @brief sets the knot_ vector and knot_values_ vector from the input map
   * @param input
   */
  void setInputData(const DataMap<T>& input);

  /**
   * @brief calculate the weight vector
   */
  void calculateWeights();

  /**
   * @brief calculates the interpolated values for the input matrix of input vectors
   * @param eval_pts: the vectors at which to interpolate
   * @return vector of interpolated values corresponding to the input vectors
   */
  std::vector<T> calculateOutput(const Eigen::Ref<RBFMatrixX<T>>& eval_pts) const;

protected:

  /** @brief radial basis function */
  typename RBFBase<T>::Ptr rbf_;

  /** @brief n_ the number of input RBF datapoints */
  std::size_t n_;

  /** @brief d_ the dimensionality of the input knots */
  std::size_t d_;

  /** @brief knots_ an n x dimensionality matrix of the input RBF knots */
  RBFMatrixX<T> knots_;

  /** @brief knot_values_ a column vector of the scores associated with each corresponding knot */
  RBFMatrixX<T> knot_values_;

  RBFMatrixX<T> rbf_matrix_;

  RBFVectorX<T> weights_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_SOLVER_H
