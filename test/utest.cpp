#include "rbf_interpolation/rbf_solver.h"
#include "rbf_interpolation/impl/rbf/gaussian_rbf.hpp"
#include "rbf_interpolation/impl/rbf/inverse_multiquadric_rbf.hpp"
#include "rbf_interpolation/impl/rbf/multiquadric_rbf.hpp"
#include "rbf_interpolation/impl/rbf/polyharmonic_rbf.hpp"
#include "rbf_interpolation/impl/rbf/thin_plate_spline_rbf.hpp"

#include <gtest/gtest.h>
#include <ros/console.h>

namespace rbf_interpolation
{
using T = double;

T surface(T x, T y)
{
  return std::pow(x, 2.0) + std::pow(y, 2.0);
}

DataMap<T> generateKnots(const T range)
{
  std::size_t n = 10;

  DataMap<T> out;
  out.reserve(n*n);

  for(std::size_t i = 0; i < n; ++i)
  {
    T x = -range / 2.0 + static_cast<T>(i / n) * range;
    for(std::size_t j = 0; j < n; ++j)
    {
      T y = -range / 2.0 + static_cast<T>(j / n) * range;
      T z = surface(x, y);
      RBFVectorX<T> knot;
      knot.resize(2);
      knot << x, y;
      out.push_back(std::pair<T, RBFVectorX<T>>(z, knot));
    }
  }

  return out;
}

RBFMatrixX<T> generateEvalPoints(T range)
{
  std::size_t n = 100;

  RBFMatrixX<T> out (n*n, 2);

  for(std::size_t i = 0; i < n; ++i)
  {
    T x = -range / 2.0 + static_cast<T>(i / n) * range;
    for(std::size_t j = 0; j < n; ++j)
    {
      T y = -range / 2.0 + static_cast<T>(j / n) * range;
      out.row(i * n + j) << x, y;
    }
  }

  return out;
}

RBFMatrixX<T> mapToMatrix(const DataMap<T>& map)
{
  RBFMatrixX<T> matrix (map.size(), map.front().second.size());

  for(std::size_t i = 0; i < map.size(); ++i)
  {
    const std::pair<const T, RBFVectorX<T>>& pair = map[i];
    matrix.row(i) = pair.second.transpose();
  }

  return matrix;
}

void testRBF(const RBFBase<T>::Ptr& rbf)
{
  RBFSolver<T> rbf_solver(rbf);

  // Sample some knots on the test surface
  T range = static_cast<T>(10.0);
  DataMap<T> knots = generateKnots(range);
  rbf_solver.setInputData(knots);

  // Calculate
  rbf_solver.calculateWeights();

  // Create points at which to evaluate the RBF representation of the surface
  RBFMatrixX<T> eval_pts = generateEvalPoints(range);
  std::vector<T> out = rbf_solver.calculateOutput(eval_pts);

  ASSERT_EQ(out.size(), eval_pts.rows()) << "Output is not the same size as the number of input evaluation points";

  // Check that the error at the knots is zero
  RBFMatrixX<T> knots_matrix = mapToMatrix(knots);
  std::vector<T> test = rbf_solver.calculateOutput(knots_matrix);

  std::size_t i = 0;
  for(const std::pair<T, RBFVectorX<T>>& pair : knots)
  {
    T error = pair.first - test[i];
    ASSERT_LT(error, std::numeric_limits<T>::epsilon()) << "Error at knot " << i << " is " << error;
    ++i;
  }
}

TEST(RBFTests, gaussianRBF)
{
  GaussianRBF<T>::Ptr rbf = GaussianRBF<T>::Ptr(new GaussianRBF<T>(1.0));
  testRBF(rbf);
}

TEST(RBFTests, inverseMultiquadricRBF)
{
  InverseMultiquadricRBF<T>::Ptr rbf = InverseMultiquadricRBF<T>::Ptr(new InverseMultiquadricRBF<T>(1.0, 1.0));
  testRBF(rbf);
}

TEST(RBFTests, multiquadricRBF)
{
  MultiquadricRBF<T>::Ptr rbf = MultiquadricRBF<T>::Ptr(new MultiquadricRBF<T>(1.0, 1.0));
  testRBF(rbf);
}

TEST(RBFTests, polyharmonicRBF)
{
  PolyharmonicRBF<T>::Ptr rbf = PolyharmonicRBF<T>::Ptr(new PolyharmonicRBF<T>(1.0));
  testRBF(rbf);
}

TEST(RBFTests, thinPlateSplineRBF)
{
  ThinPlateSplineRBF<T>::Ptr rbf = ThinPlateSplineRBF<T>::Ptr(new ThinPlateSplineRBF<T>());
  testRBF(rbf);
}

} // namespace rbf_interpolation

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  ros::Time::init();
  return RUN_ALL_TESTS();
}
