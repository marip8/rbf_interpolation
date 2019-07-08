#ifndef RBF_INTERPOLATION_POLYHARMONIC_RBF_H
#define RBF_INTERPOLATION_POLYHARMONIC_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class phi(x) = (||x-c||)^k
 */
template<typename T>
class PolyharmonicRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<PolyharmonicRBF> Ptr;
  using Vector = typename RBFVectorX<T>::type;

  PolyharmonicRBF(const T k)
    : k_(k)
  {

  }

  T calculate(const Vector& pt, const Vector& center) const override
  {
    Vector vec = pt - center;
    T r = vec.norm();
    return std::pow(r, k_);
  }

protected:

  T k_;
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_POLYHARMONIC_RBF_H
