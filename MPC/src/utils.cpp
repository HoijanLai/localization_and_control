#include "utils.h"
#include <math.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include <cppad/cppad.hpp>
using CppAD::AD;


// Evaluate a polynomial.
AD<double> polyeval(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}

// Evaluate a polynomial derivative
AD<double> polyderivative(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * CppAD::pow(x, i - 1);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

// transform form global to local coordinate
Eigen::MatrixXd global2local(std::vector<double> ptsx, std::vector<double> ptsy, double psi, double px, double py) {
  int length = ptsx.size();
  Eigen::MatrixXd result(2, length);
  for (int i = 0; i < length; ++i) {
    double x = ptsx[i] - px;
    double y = ptsy[i] - py;
    result(0, i) =   x * cos(psi) + y * sin(psi);
    result(1, i) = - x * sin(psi) + y * cos(psi);
  }
  return result;
}

// Eigen library to std
std::vector<double> eigen2std(Eigen::VectorXd eigen_vec) {
  std::vector<double> result;
  result.resize(eigen_vec.size());
  Eigen::VectorXd::Map(&result[0], eigen_vec.size()) = eigen_vec;
  return result;
}
