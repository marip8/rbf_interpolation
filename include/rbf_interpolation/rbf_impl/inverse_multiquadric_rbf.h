#ifndef RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H
#define RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the form \f$\phi(x) = 1 / (1 + (\epsilon * ||x-c||)^2)^(k/2)\f$
 */
template<typename T>
class InverseMultiquadricRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<InverseMultiquadricRBF<T>> Ptr;

  inline explicit InverseMultiquadricRBF(const T eps = static_cast<T>(1.0),
                                         const T k = static_cast<T>(1.0))
    : eps_(std::abs(eps))
    , k_(std::abs(k))
  {
  }

  inline virtual T calculate(const Eigen::Ref<RBFVectorX<T>>& pt,
                             const Eigen::Ref<RBFVectorX<T>>& center) const override
  {
    RBFVectorX<T> vec = pt - center;
    T r = vec.norm();
    T val = static_cast<T>(1) + std::pow(eps_ * r, static_cast<T>(2.0));
    T denom = std::pow(val, k_ / static_cast<T>(2.0));
    return (static_cast<T>(1.0) / denom);
  }

  inline virtual unsigned order() const override
  {
    return 0;
  }

private:

  /** @brief shape parameter on (0, inf) */
  T eps_;
  /** @brief order on [1, inf) */
  T k_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H
