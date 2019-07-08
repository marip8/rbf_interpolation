#ifndef RBF_INTERPOLATION_RBF_BASE_H
#define RBF_INTERPOLATION_RBF_BASE_H

#include <memory>
#include "rbf_matrix_types.h"

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
  using Vector = typename RBFVectorX<T>::type;

  RBFBase() = default;

  virtual ~RBFBase()
  {

  }

  virtual T calculate(const Vector& pt, const Vector& center) const = 0;

};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_RBF_BASE_H
