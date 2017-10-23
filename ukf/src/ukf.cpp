#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


#define EPS 0.00001

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initialization flag
  is_initialized_ = false;
  time_us_ = 0;

  // sensor usage triggers
  use_laser_ = true;
  use_radar_ = true;

  // Process noise
  std_a_ = 1;
  std_yawdd_ = 4;

  // Laser noise
  std_laspx_ = 0.15;
  std_laspy_ = 0.15;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << pow(std_laspx_, 2),                 0,
                               0, pow(std_laspy_, 2);

  // Radar noise
  std_radr_ = 0.3;
  std_radphi_ = 0.03;
  std_radrd_ = 0.3;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_, 2),                   0,                  0,
                              0, pow(std_radphi_, 2),                  0,
                              0,                   0, pow(std_radrd_, 2);

  // dimensions
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  // placeholders
  x_ = VectorXd::Zero(n_x_);
  P_ = MatrixXd::Identity(n_x_, n_x_);
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);

  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i)
      weights_(i) = 0.5 / (lambda_ + n_aug_);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (use_radar_ == false && use_laser_ == false) return;

  if (!is_initialized_) {
    // directly take the measurement as the current state
    time_us_   = meas_package.timestamp_;
    VectorXd z = meas_package.raw_measurements_;
    x_ << 0, 0, 8, 0, 0;
    if      (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      x_(0) = z(0) * cos(z(1));
      x_(1) = z(0) * sin(z(1));
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
      x_(0) = z(0);
      x_(1) = z(1);
    }

    is_initialized_ = true;
    return;
  }


  // get delta_t
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  //*************************
  // predict
  //*************************
  Prediction(dt);

  //************************
  // update
  //************************

  if      (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    cout << "radar" << endl;
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    cout << "laser" << endl;
    UpdateLidar(meas_package);
  }

  cout << "x: " << x_ << endl;
  cout << "P: " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {

  //*********************************************************
  //  Generate the augmented sigma points
  //*********************************************************

  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create augmented sigma point matrix
  MatrixXd Xsig_aug(n_aug_, 2*n_aug_+1);
  MatrixXd A_aug = P_aug.llt().matrixL();
  Xsig_aug.col(0) = x_aug;
  for (int i = 1; i <= n_aug_; ++i) {
      Xsig_aug.col(i)        = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(i-1);
      Xsig_aug.col(i+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(i-1);
  }


  //*********************************************************
  //  Predict the augmented sigma points
  //*********************************************************

  // modify by column
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // unpack the sigma point vector
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // initialize the delta state components
    double dpx = 0.0, dpy = 0.0, dv = 0.0, dyaw = 0.0, dyawd = 0.0;

    // transition part of the delta state components
    if (fabs(yawd) > EPS) {
        dpx += v/yawd * ( sin(yaw + yawd * dt) - sin(yaw));
        dpy += v/yawd * (-cos(yaw + yawd * dt) + cos(yaw));
    }
    else {
        dpx += v * cos(yaw) * dt;
        dpy += v * sin(yaw) * dt;
    }

    dyaw += yawd * dt;

    // noise part of the delta state components
    dpx   += 0.5 * dt * dt * cos(yaw) * nu_a;
    dpy   += 0.5 * dt * dt * sin(yaw) * nu_a;
    dv    +=            dt * nu_a;
    dyaw  += 0.5 * dt * dt * nu_yawdd;
    dyawd +=            dt * nu_yawdd;

    // add delta to the current sigma state
    Xsig_pred_(0,i) = px   + dpx;
    Xsig_pred_(1,i) = py   + dpy;
    Xsig_pred_(2,i) = v    + dv;
    Xsig_pred_(3,i) = yaw  + dyaw;
    Xsig_pred_(4,i) = yawd + dyawd;
  }


  //****************************************************************
  // get predicted state and covariance from the sigma predictions
  //****************************************************************
  // predict state mean
  x_ = Xsig_pred_*weights_;

  // predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = shrink(x_diff(3));
      P_ += weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  char n_z = 2;

  //************************************************
  //  Predict measurement
  //***********************************************

  // transform sigma points into measurement space
  MatrixXd Zsig = Xsig_pred_.topRows(n_z);
  VectorXd z_pred = Zsig * weights_;

  // calculate measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R_laser_;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff =       Zsig.col(i) - z_pred;
      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // Calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS:" << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  char n_z = 3;

  //****************************************************
  //  Predict measurement
  //****************************************************

  // transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2*n_aug_+1);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      // unpack the sigma prediction vector
      double px   = Xsig_pred_(0, i);
      double py   = Xsig_pred_(1, i);
      double v    = Xsig_pred_(2, i);
      double yaw  = Xsig_pred_(3, i);

      double rho = sqrt(px*px + py*py);
      if (fabs(rho) > EPS)
          Zsig.col(i) << rho,
                         atan2(py, px),
                         (px*cos(yaw)*v + py*sin(yaw)*v) / rho;
  }

  // calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;


  //****************************************************
  //  Update
  //****************************************************

  // calculate measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      VectorXd z_diff = Zsig.col(i) - z_pred;

      z_diff(1) = shrink(z_diff(1));
      S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R_radar_;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff =       Zsig.col(i) - z_pred;

      x_diff(3) = shrink(x_diff(3));
      z_diff(1) = shrink(z_diff(1));

      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  z_diff(1) = shrink(z_diff(1));

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS:" << NIS_radar_ << endl;

}



/**
 * normalized angle, so it between -pi and pi
 * @param the angle to normalized
 */
double UKF::shrink(double angle) {
  double ang = angle;
  if (angle > M_PI) {
    double r = fmod((angle - M_PI), (2 * M_PI));
    ang = r - M_PI;
  }
  if (angle < -M_PI) {
    double r = fmod((angle + M_PI) ,(2 * M_PI));
    ang = r + M_PI;
  }
  return ang;
}
