#ifndef RBF_INTERPOLATION_POLYHARMONIC_RBF_H
#define RBF_INTERPOLATION_POLYHARMONIC_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class \f$\phi(x) = (||x-c||)^(2*k - 1)\f$
 */
template<typename T>
class PolyharmonicRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<PolyharmonicRBF> Ptr;

  inline explicit PolyharmonicRBF(const T k = static_cast<T>(1.0))
    : k_(std::abs(k))
  {

  }

  inline virtual T calculate(const Eigen::Ref<RBFVectorX<T>>& pt,
                             const Eigen::Ref<RBFVectorX<T>>& center) const override
  {
    RBFVectorX<T> vec = pt - center;
    T r = vec.norm();
    T exponent = static_cast<T>(2.0) * k_ - static_cast<T>(1.0);
    return std::pow(r, exponent);
  }

  inline virtual unsigned order() const override
  {
    return std::ceil(k_);
  }

protected:
  /** @brief order on [1, inf) */
  T k_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_POLYHARMONIC_RBF_H
