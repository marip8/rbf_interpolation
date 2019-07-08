#ifndef RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H
#define RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class phi(x) = (||x-c||)^2*ln(||x-c||)
 */
template<typename T>
class ThinPlateSplineRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<ThinPlateSplineRBF<T>> Ptr;
  using Vector = typename RBFVectorX<T>::type;

  ThinPlateSplineRBF()
  {

  }

  T calculate(const Vector& pt, const Vector& center) const override
  {
    Vector vec = pt - center;
    T r = vec.norm();
    return std::pow(r, static_cast<T>(2.0)) * std::log(r);
  }
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H
