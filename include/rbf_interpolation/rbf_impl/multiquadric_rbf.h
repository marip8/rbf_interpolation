#ifndef RBF_INTERPOLATION_MULTIQUADRIC_RBF_H
#define RBF_INTERPOLATION_MULTIQUADRIC_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the form phi = sqrt(1 + (epsilon * ||x-c||)^2)
 */
template<typename T>
class MultiquadricRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<MultiquadricRBF<T>> Ptr;
  using Vector = typename RBFVectorX<T>::type;

  MultiquadricRBF(const T eps)
    : eps_(eps)
  {

  }

  T calculate(const Vector& pt, const Vector& center) const override
  {
    Vector vec = pt - center;
    T r = vec.norm();
    return std::sqrt(static_cast<T>(1) + std::pow(eps_ * r, static_cast<T>(2.0)));
  }

private:

  T eps_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_MULTIQUADRIC_RBF_H
