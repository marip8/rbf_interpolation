#ifndef RBF_INTERPOLATION_RBF_BASE_H
#define RBF_INTERPOLATION_RBF_BASE_H

#include <memory>
#include "rbf_interpolation/rbf_matrix_types.h"

namespace rbf_interpolation
{

/**
 * @brief Base class for radial basis functions
 */
template<typename T>
class RBFBase
{
public:
  typedef typename std::shared_ptr<RBFBase<T>> Ptr;

  RBFBase() = default;
  virtual ~RBFBase() = default;

  /**
   * @brief evaluates the radial basis function given two input points
   * @param pt
   * @param center
   * @return
   */
  virtual T calculate(const Eigen::Ref<RBFVectorX<T>>& pt,
                      const Eigen::Ref<RBFVectorX<T>>& center) const = 0;

  /**
   * @brief returns the order of the required polynomial tail for the radial basis function
   */
  virtual unsigned order() const = 0;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_BASE_H
