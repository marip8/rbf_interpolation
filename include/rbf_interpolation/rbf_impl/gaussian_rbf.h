#ifndef RBF_INTERPOlATION_GAUSSIAN_RBF_H
#define RBF_INTERPOlATION_GAUSSIAN_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class phi(x) = e^(-(epslion * ||x-c||)^2)
 * @param r0
 */
template<typename T>
class GaussianRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<GaussianRBF<T>> Ptr;
  using Vector = typename RBFVectorX<T>::type;

  GaussianRBF(const T eps)
    : eps_(eps)
  {

  }

  T calculate(const Vector& pt, const Vector& center) const override
  {
    Vector vec = pt - center;
    T r = vec.norm();
    return std::exp(-std::pow(eps_ * r , static_cast<T>(2.0)));
  }

protected:

  T eps_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOlATION_GAUSSIAN_RBF_H
