#ifndef RBF_INTERPOLATION_RBF_SOLVER_IMPL_H
#define RBF_INTERPOLATION_RBF_SOLVER_IMPL_H

#include "rbf_solver.h"

namespace rbf_interpolation
{

template<typename T>
RBFSolver<T>::RBFSolver(const RBFBasePtr rbf)
  : rbf_(rbf)
{

}

template<typename T>
void RBFSolver<T>::setInputData(const RBFDataMap& input)
{
  ROS_DEBUG("Starting to save input data");

  // Resize the containers
  n_ = input.size();
  d_ = input.begin()->second.size();
  ROS_DEBUG_STREAM("Size and dimensionality: " << n_ << " x " << d_);

  knots_.resize(n_, d_);
  knot_values_.resize(n_ + d_ + 1, 1);
  std::size_t idx = 0;

  for(const auto& it : input)
  {    
    knot_values_(idx) = it.first;
    knots_.row(idx) = it.second;
    ++idx;
  }

  // Set the last d_ + 1 elements to zero
//  knot_values_.segment(n_, d_ + 1) = RBFVector::Zero(d_ + 1);
  for(std::size_t i = n_; i < n_ + d_ + 1; ++i)
  {
    knot_values_(i) = static_cast<T>(0.0);
  }

  ROS_DEBUG_STREAM("Finished saving input data"); //: \n" << knot_values_);
}

template<typename T>
void RBFSolver<T>::calculateWeights()
{
  ROS_DEBUG("Started calculating weights");

  ros::Time start = ros::Time::now();

  /* Create the RBF matrix to solve the following equation
   * | K      P | * | w | = | k |
   * | P^T    0 |   |   |   | 0 |
   *
   * where:
   * K(i,j) = rbf(knot(i), knot(j))
   * P(i) = [1, x1, x2, ..., xn]
   * w is the vector of weights to solve for
   * [k, 0]^T are the knot values and vector of zeros of size (d_ + 1)
   */

  // Resize the matrix
  const std::size_t size = n_ + d_ + 1;
  rbf_matrix_.resize(size, size);
  weights_.resize(size);

  // Loop over the input points and evaluate the RBF for each combination of data inputs
  for(std::size_t i = 0; i < n_; ++i)
  {
    for(std::size_t j = 0; j < n_; ++j)
    {
      rbf_matrix_(i, j) = rbf_->calculate(knots_.row(i), knots_.row(j));
    }
  }

  // Drop in the one vectors and input locators for the polynomial tail
  rbf_matrix_.block(0, n_, n_, 1) = RBFVector::Ones(n_);
  rbf_matrix_.block(n_, 0, 1, n_) = RBFVector::Ones(n_).transpose();

  rbf_matrix_.block(0, n_ + 1, n_, d_) = knots_;
  rbf_matrix_.block(n_ + 1, 0, d_, n_) = knots_.transpose();

  // Set the bottom right (d_ + 1) x (d_ + 1) block to all zeros
  rbf_matrix_.bottomRightCorner(d_ + 1, d_ + 1).setZero();

  // Get the weights by solving the matrix equation with the given input data
  // The RBF matrix should be positive definite, so a solver like LLT or LLDT should work
//    weights_ = rbf_matrix_.llt().solve(knot_values_);
//  weights_ = rbf_matrix_.ldlt().solve(knot_values_);
  weights_ = rbf_matrix_.fullPivLu().solve(knot_values_);

  ROS_DEBUG_STREAM("Finished calculating weights (" << (ros::Time::now() - start).toSec()<< " sec)"); //: \n" << weights_);
}

template<typename T>
std::vector<T> RBFSolver<T>::calculateOutput(const RBFMatrix& eval_pts) const
{
  ROS_DEBUG("Starting output calculation");

  // Initialize a results vector
  std::vector<T> output;

  if(weights_.size() == 0)
  {
    ROS_ERROR("No weights have been calculated yet");
    return {};
  }
  else
  {
    for(std::size_t i = 0; i < eval_pts.rows(); ++i)
    {
      const RBFVector& pt = eval_pts.row(i);

      // Create a column vector for the evaluation of the RBF at each input point
      RBFVector r;
      r.resize(n_ + d_ + 1);

      for(std::size_t j = 0; j < n_; ++j)
      {
        const RBFVector& center = knots_.row(j);

        if(pt.size() != center.size())
        {
          ROS_ERROR("Point and center are not the same size");
          return {};
        }

        r[j] = rbf_->calculate(pt, center);
      }

      r[n_] = 1.0;
      r.segment(n_+ 1, d_) = pt.transpose();

      // Multiply the weight vector by the r vector to get the value of the function at the given input point
      T result = weights_.transpose() * r;
      output.push_back(std::move(result));
    }
  }

  ROS_DEBUG("Finished output calculation");

  return output;
}

} // namespace rbf_interpolation


#endif // RBF_INTERPOLATION_RBF_SOLVER_IMPL_H
