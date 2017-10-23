#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // size check
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    cout << "Invalid estimations or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
      VectorXd residuals = estimations[i] - ground_truth[i];
      residuals = residuals.array() * residuals.array();
      rmse += residuals;
  }

  // take the average
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;

}
