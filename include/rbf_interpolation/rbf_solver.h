#ifndef RBF_INTERPOLATION_RBF_SOLVER_H
#define RBF_INTERPOLATION_RBF_SOLVER_H

#include <map>
#include <functional>
#include <memory>
#include "rbf_base.h"
#include <ros/console.h>

namespace rbf_interpolation
{

template <typename T>
class RBFSolver
{
public:
  using RBFVector = typename RBFVectorX<T>::type;
  using RBFMatrix = typename RBFMatrixX<T>::type;
  using RBFDataMap = typename DataMap<T>::type;
  typedef typename RBFBase<T>::Ptr RBFBasePtr;

  RBFSolver(const RBFBasePtr rbf);

  void setInputData(const RBFDataMap& input);

  void calculateWeights();

  std::vector<T> calculateOutput(const RBFMatrix& eval_pts) const;

protected:

  RBFBasePtr rbf_;

  /**
   * @brief n_ the number of input RBF datapoints
   */
  std::size_t n_;

  /**
   * @brief d_ the dimensionality of the input knots
   */
  std::size_t d_;

  /**
   * @brief knots_ an n x dimensionality matrix of the input RBF knots
   */
  RBFMatrix knots_;

  /**
   * @brief knot_values_ a column vector of the scores associated with each corresponding knot
   */
  RBFMatrix knot_values_;

  RBFMatrix rbf_matrix_;

  RBFVector weights_;

};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_SOLVER_H
