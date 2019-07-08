#include <rbf_interpolation/rbf_solver_impl.h>
#include <rbf_interpolation/rbf_impl/gaussian_rbf.h>
#include <rbf_interpolation/rbf_impl/inverse_multiquadric_rbf.h>
#include <rbf_interpolation/rbf_impl/multiquadric_rbf.h>
#include <rbf_interpolation/rbf_impl/polyharmonic_rbf.h>
#include <rbf_interpolation/rbf_impl/thin_plate_spline_rbf.h>

#include <gtest/gtest.h>
#include <ros/ros.h>

namespace rbf_interpolation
{
using T = double;
using RBFDataMap = DataMap<T>::type;
using RBFVector = RBFVectorX<T>::type;
using RBFMatrix = RBFMatrixX<T>::type;

T surface(T x, T y)
{
  return std::pow(x, 2.0) + std::pow(y, 2.0);
}

RBFDataMap generateKnots(const T range)
{
  std::size_t n = 10;

  RBFDataMap out;
  out.reserve(n*n);

  for(std::size_t i = 0; i < n; ++i)
  {
    T x = -range / 2.0 + static_cast<T>(i / n) * range;
    for(std::size_t j = 0; j < n; ++j)
    {
      T y = -range / 2.0 + static_cast<T>(j / n) * range;
      T z = surface(x, y);
      RBFVector knot;
      knot.resize(2);
      knot << x, y;
      out.push_back(std::pair<T, RBFVector>(z, knot));
    }
  }

  return out;
}

RBFMatrix generateEvalPoints(T range)
{
  std::size_t n = 100;

  RBFMatrix out (n*n, 2);

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

RBFMatrix mapToMatrix(const RBFDataMap& map)
{
  RBFMatrix matrix (map.size(), map.front().second.size());

  for(std::size_t i = 0; i < map.size(); ++i)
  {
    const std::pair<const T, RBFVector>& pair = map[i];
    matrix.row(i) = pair.second.transpose();
  }

  return matrix;
}

void testRBF(const RBFBase<T>::Ptr& rbf)
{
  RBFSolver<T> rbf_solver(rbf);

  // Sample some knots on the test surface
  T range = static_cast<T>(10.0);
  RBFDataMap knots = generateKnots(range);
  rbf_solver.setInputData(knots);

  // Calculate
  rbf_solver.calculateWeights();

  // Create points at which to evaluate the RBF representation of the surface
  RBFMatrix eval_pts = generateEvalPoints(range);
  std::vector<T> out = rbf_solver.calculateOutput(eval_pts);

  ASSERT_EQ(out.size(), eval_pts.rows()) << "Output is not the same size as the number of input evaluation points";

  // Check that the error at the knots is zero
  RBFMatrix knots_matrix = mapToMatrix(knots);
  std::vector<T> test = rbf_solver.calculateOutput(knots_matrix);

  std::size_t i = 0;
  for(const std::pair<T, RBFVector>& pair : knots)
  {
    T error = pair.first - test[i];
    ASSERT_TRUE(error < static_cast<T>(1e-6)) << "Error at knot " << i << " is " << error;
    ++i;
  }
}

TEST(RBFTests, gaussianRBF)
{
  GaussianRBF<T>::Ptr rbf = GaussianRBF<T>::Ptr(new GaussianRBF<T>(10.0));
  testRBF(rbf);
}

TEST(RBFTests, inverseMultiquadricRBF)
{
  InverseMultiquadricRBF<T>::Ptr rbf = InverseMultiquadricRBF<T>::Ptr(new InverseMultiquadricRBF<T>(10.0));
  testRBF(rbf);
}

TEST(RBFTests, multiquadricRBF)
{
  MultiquadricRBF<T>::Ptr rbf = MultiquadricRBF<T>::Ptr(new MultiquadricRBF<T>(10.0));
  testRBF(rbf);
}

TEST(RBFTests, polyharmonicRBF)
{
  PolyharmonicRBF<T>::Ptr rbf = PolyharmonicRBF<T>::Ptr(new PolyharmonicRBF<T>(3.0));
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
  ros::init(argc, argv, "rbf_interpolation_test");
  ros::Time::init();
  return RUN_ALL_TESTS();
}
