#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
# define PI  3.14159265358979323846  /* pi */
# define EPS 1e-5

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;

  _update_state_cov(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  double rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  double phi = atan2(x_(1), x_(0));
  double rho_d = 0;

  if (fabs(rho) >= EPS)
    rho_d = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_d;
  VectorXd y = z - z_pred;

  // map the phi between [-pi, pi]
  if (y(1) >  PI) y(1) -= 2*PI;
  if (y(1) < -PI) y(1) += 2*PI;

  _update_state_cov(y);
}


void KalmanFilter::_update_state_cov(const VectorXd &y) {
  /**
    * Common calculations in updating
  */
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
