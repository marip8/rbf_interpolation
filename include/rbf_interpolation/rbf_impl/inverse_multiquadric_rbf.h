#ifndef RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H
#define RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the form phi = 1 / sqrt(1 + (epsilon * ||x-c||)^2)
 */
template<typename T>
class InverseMultiquadricRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<InverseMultiquadricRBF<T>> Ptr;
  using Vector = typename RBFVectorX<T>::type;

  InverseMultiquadricRBF(const T eps)
    : eps_(eps)
  {

  }

  T calculate(const Vector& pt, const Vector& center) const override
  {
    Vector vec = pt - center;
    T r = vec.norm();
    return (static_cast<T>(1.0) / std::sqrt(static_cast<T>(1) + std::pow(eps_ * r, static_cast<T>(2.0))));
  }

private:

  T eps_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_INVERSE_MULTIQUADRIC_RBF_H
