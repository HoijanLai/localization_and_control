#ifndef UTILS_H
#define UTILS_H

#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>
using CppAD::AD;

AD<double> polyeval(Eigen::VectorXd coeffs, AD<double> x);

AD<double> polyderivative(Eigen::VectorXd coeffs, AD<double> x);

Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);

Eigen::MatrixXd global2local(std::vector<double> ptsx, std::vector<double> ptsy, double psi, double px, double py);

std::vector<double> eigen2std(Eigen::VectorXd eigen_vec);

#endif
