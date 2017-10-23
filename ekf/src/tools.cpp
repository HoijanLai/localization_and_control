#include <iostream>
#include "tools.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

#define EPS 1e-5


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
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

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /* calculate the Jacobian */
    MatrixXd Hj(3, 4);

    // recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    // denominators
    double c1 = px*px + py*py;
    double c2 = sqrt(c1);
    double c3 = (c1*c2);

    // check division by zero
    if (fabs(c1) < EPS) {
      std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
      return Hj;
    }

    Hj <<            (px / c2),               (py / c2),       0,       0,
                    -(py / c1),               (px / c1),       0,       0,
       py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

   return Hj;
}
