#include "kalman_filter.h"
#include <iostream>
#include "math.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * Predict the state
   */
   x_ = F_ * x_;
   P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update the state by using Kalman Filter equations
   */
   VectorXd z_pred = H_ * x_;
   VectorXd y = z - z_pred;
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  // Compute the h'(x) matrix coefficients
  float rho_ = sqrt(x_(0)*x_(0) + x_(1)*x_(1));

  // avoid dividing by zero. Arbitrarily, if rho is too small, we set it
  // to a very small value. This prevent dividing by zero
  if (rho_ < 0.00001) {
    rho_ = 0.00001;
  }
  // rest of h'(x) coefficients
  float phi_ = atan2(x_(1),x_(0));
  float rhodot_ = (x_(0)*x_(2) + x_(1)*x_(3)) / rho_;

  // Compute y = z_radar - h'(x)
  // h'(x) is called z_pred here
  VectorXd z_pred(3);
  z_pred << rho_, phi_, rhodot_;
  VectorXd y = z - z_pred;

  // We need to normalize angle as change from -pi to pi can be problematic:
  // If the angle is greater than Pi, we can substract 2*pi until the angle
  // is in the range [-pi; pi]
  while (y(1) > M_PI) {
    y(1) -= (2 * M_PI);
  }
  // Likewise, if the angle is smaller than -pi, we add 2*pi until the angle
  // is in the range [-pi:pi]
  while (y(1) < -M_PI) {
    y(1) += (2 * M_PI);
  }


  // Back to the same steps than LiDar update
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;


  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
