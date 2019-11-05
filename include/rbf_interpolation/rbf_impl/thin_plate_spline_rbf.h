#ifndef RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H
#define RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H

#include <rbf_interpolation/rbf_base.h>

namespace rbf_interpolation
{

/**
 * @brief Radial basis function of the class \f$\phi(x) = (||x-c||)^2*ln(||x-c||)\f$
 */
template<typename T>
class ThinPlateSplineRBF : public RBFBase<T>
{
public:
  typedef typename std::shared_ptr<ThinPlateSplineRBF<T>> Ptr;

  inline explicit ThinPlateSplineRBF() = default;

  inline virtual T calculate(const Eigen::Ref<RBFVectorX<T>>& pt,
                             const Eigen::Ref<RBFVectorX<T>>& center) const override
  {
    RBFVectorX<T> vec = pt - center;
    const T r = vec.norm();
    const T result = std::pow(r, static_cast<T>(2.0)) * std::log(r);

    return std::isnan(result) ? static_cast<T>(0.0) : result;
  }

  inline virtual unsigned order() const override
  {
    return 2;
  }
};

} // namespace rbf_interpolation

#endif // RBF_INTERPOLATION_THIN_PLATE_SPLINE_RBF_H
