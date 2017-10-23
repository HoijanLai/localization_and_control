#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);


  // H_laser_remains constant
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225,      0,
                   0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  // placeholder for [F], [P] and [H_j]
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0,    0,    0,
             0, 1,    0,    0,
             0, 0, 1000,    0,
             0, 0,    0, 1000;

  Hj_ = MatrixXd(3, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}




void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
   /**
   the first generation of state x is just the measurement.
   */
   if (!is_initialized_) {


     // first measurement
     cout << "EKF: " << endl;
     ekf_.x_ = VectorXd(4);
     ekf_.x_ << 1, 1, 1, 1;

     previous_timestamp_ = measurement_pack.timestamp_;

     if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

       // cout << "radar" << endl;
       Eigen::VectorXd z = measurement_pack.raw_measurements_;
       // polar to cartesian
       ekf_.x_ << z(0)*cos(z(1)),
                  z(0)*sin(z(1)),
                  z(2)*cos(z(1)),
                  z(2)*sin(z(1));
     }

     else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
       // cout << "laser" << endl;
       Eigen::VectorXd z = measurement_pack.raw_measurements_;
       ekf_.x_ << z(0), z(1), 0, 0;
     }

       // done initializing, no need to predict or update
       // cout << "initialized" << endl;
       is_initialized_ = true;
       return;
   }




  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

   double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

   ekf_.F_(0, 2) = dt;
   ekf_.F_(1, 3) = dt;


   double noise_ax = 9;
   double noise_ay = 9;
   double dt2 = dt*dt;
   double dt3 = dt2*dt;
   double dt4 = dt3*dt;

   ekf_.Q_ = MatrixXd(4, 4);
   ekf_.Q_ << dt4*noise_ax/4,              0, dt3*noise_ax/2,              0,
                           0, dt4*noise_ay/4,              0, dt3*noise_ay/2,
              dt3*noise_ax/2,              0,   dt2*noise_ax,              0,
                           0, dt3*noise_ay/2,              0,   dt2*noise_ay;


  ekf_.Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    printf("got radar measurement\n");
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }

  else {
    // Laser updates
    printf("got laser measurement\n");
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
