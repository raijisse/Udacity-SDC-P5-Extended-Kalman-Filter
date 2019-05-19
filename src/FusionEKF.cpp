#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
   // Laser projection matrix
   H_laser_ << 1, 0, 0, 0,
               0, 1, 0, 0;

   // State covariance matrix P
   MatrixXd P_ = MatrixXd(4, 4);
   P_ << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

   // The initial transition matrix F_
   MatrixXd F_ = MatrixXd(4,4);
   F_ << 1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0,
         0, 0, 0, 1;

   // Initialize the process noise covariance matrix Q
   MatrixXd Q_ = MatrixXd(4, 4);

   // Initialize the state vector
   VectorXd x_ = VectorXd(4);

   // Finally, initialize the Kalman filter with the previously set
   // values
   ekf_.Init(x_, P_,F_,H_laser_,R_laser_,Q_);

}


/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     */

    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // If the first data is radar data:
      // We initialize the phi, rho, and rhodot measurements
      float phi_ = measurement_pack.raw_measurements_[0];
      float rho_ = measurement_pack.raw_measurements_[1];
      float rhodot_ = measurement_pack.raw_measurements_[2];

      // Update state with cartesian coordinates derived from radar data
      ekf_.x_ << cos(phi_),
                 sin(phi_),
                 0,
                 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // If first data is LiDar measurements:
      // We initialize state with raw lidar measurements

      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
    }

    // Store timestamp from first observation for computing dt later on
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Update the process noise covariance matrix.
   */

  // set the acceleration noise components
  int noise_ax = 9;
  int noise_ay = 9;

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Compute all dt derivative for later computations
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the process noise covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // Predict
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // If data is from radar

    // Reset the R_ for the appropriate sensor for computations
    ekf_.R_ = R_radar_;
    // Compute Jacobian from sensor state
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    // Update using Extended KF and radar data
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);


  } else {
    // If data is from LiDar

    // Reset the R_ and H_ for the appropriate sensor for computations
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    // Update using KF and LiDar data
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
