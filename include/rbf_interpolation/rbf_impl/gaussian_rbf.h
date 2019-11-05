#ifndef RBF_INTERPOlATION_GAUSSIAN_RBF_H
#define RBF_INTERPOlATION_GAUSSIAN_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class \f$\phi(x) = e^(-(\epsilon * ||x-c||)^2)\f$
 */
template<typename T>
class GaussianRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<GaussianRBF<T>> Ptr;

  inline explicit GaussianRBF(const T eps = static_cast<T>(1.0))
    : eps_(std::abs(eps))
  {
  }

  inline virtual T calculate(const Eigen::Ref<RBFVectorX<T>>& pt,
                             const Eigen::Ref<RBFVectorX<T>>& center) const override
  {
    RBFVectorX<T> vec = pt - center;
    T r = vec.norm();
    return std::exp(-std::pow(eps_ * r , static_cast<T>(2.0)));
  }

  inline virtual unsigned order() const override
  {
    return 1;
  }

protected:
  /** @brief shape parameter on (0, inf)*/
  T eps_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOlATION_GAUSSIAN_RBF_H
