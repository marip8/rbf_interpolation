#include "rbf_interpolation/impl/rbf/gaussian_rbf.hpp"
#include "rbf_interpolation/impl/rbf/inverse_multiquadric_rbf.hpp"
#include "rbf_interpolation/impl/rbf/multiquadric_rbf.hpp"
#include "rbf_interpolation/impl/rbf/polyharmonic_rbf.hpp"
#include "rbf_interpolation/impl/rbf/thin_plate_spline_rbf.hpp"

namespace rbf_interpolation
{

template class GaussianRBF<float>;
template class GaussianRBF<double>;

template class InverseMultiquadricRBF<float>;
template class InverseMultiquadricRBF<double>;

template class MultiquadricRBF<float>;
template class MultiquadricRBF<double>;

template class PolyharmonicRBF<float>;
template class PolyharmonicRBF<double>;

template class ThinPlateSplineRBF<float>;
template class ThinPlateSplineRBF<double>;

} // rbf_interpolation

