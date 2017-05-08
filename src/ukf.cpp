#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

using std::domain_error;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Previous measurement time stamp
  previous_timestamp_ = 0;

  // Keep track of initialization
  initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set measurement dimension
  n_z_ = 3;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create vector for predicted state
  x_ = VectorXd(n_x_);
  x_.fill(0);

  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);
  Xsig_.fill(0);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0);

  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig_.fill(0);

  z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0);
  
  S_ = MatrixXd(n_z_, n_z_);
  S_.fill(0);

  //create covariance matrix for prediction
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0);

  R_ = MatrixXd(n_z_, n_z_);
  R_ << std_radr_ * std_radr_, 0.0, 0.0,
        0.0, std_radphi_ * std_radphi_, 0.0,
        0.0, 0.0, std_radrd_ * std_radrd_;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  
  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  double delta_t = 0.0;
  double px, py, v, psi, psi_dot;
  VectorXd measurement_matrix;
  
  // Extract measurement data
  try {
    measurement_matrix = ExtractMeasurementData(meas_package, previous_timestamp_, delta_t);
  } catch (domain_error e) {
    // Something was wrong with the measurement package
    // Log the error and exit the function.
    cout << e.what() << endl;
    return;
  }

  // If we have not initialized our matrices
  // Just initialize and skip prediction/update until the next cycle
  if (!initialized_) {    
    px = measurement_matrix(0);
    py = measurement_matrix(1);
    v = 1.0;
    psi = 1.0;
    psi_dot = 1.0;

    x_ << px, py, v, psi, psi_dot; 
    P_.setIdentity();

    previous_timestamp_ = meas_package.timestamp_;

    initialized_ = true;
    
    return;
  }

  // Perform the prediction step
  Prediction(delta_t);

  // Perform the update step
  Update(meas_package);
}

/**
 * Extracts the appropriate measurement data and returns the measurement vector
 * @param {MeasurementPackage} meas_package
 * @param {double &} delta_t will be updated with the elapsed time after calculation
 */
VectorXd UKF::ExtractMeasurementData(MeasurementPackage meas_package, long long &previous_timestamp, double &delta_t) {
  VectorXd measurement_matrix = VectorXd(2);
  double px, py;
  double rho, rho_dot, phi;

  if (!use_radar_ && !use_laser_) {
    throw domain_error("UKF::ExtractMeasurementData: Both use_radar_ and use_laser_ are False. No data to process.");
  }

  // Update the elapsed time and the previous timestamp
  delta_t = (meas_package.timestamp_ - previous_timestamp) / 1000000.0;
  previous_timestamp = meas_package.timestamp_;

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    rho = meas_package.raw_measurements_(0);
    phi = meas_package.raw_measurements_(1);
    rho_dot = meas_package.raw_measurements_(2);

    px = rho * cos(phi);
    py = rho * sin(phi);

    measurement_matrix << px, py;

  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {

    px = meas_package.raw_measurements_(0);
    py = meas_package.raw_measurements_(1);

    measurement_matrix << px, py;

  }

  return measurement_matrix;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //create vector for weights
  double w;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0);

  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    if (i == 0) {
      w = lambda_ / (lambda_ + n_aug_);          
    } else {
      w = 1 / (2 * (lambda_ + n_aug_));
    }
    
    weights_(i) = w;
  }


  //create example matrix with predicted sigma points

  // Generate sigma points
  GenerateSigmaPoints();

  // Augment sigma points with process noise
  MatrixXd Xsig_aug = AugmentSigmaPoints();

  // Predict sigma points
  PredictSigmaPoints(Xsig_aug, delta_t);

  // Transform sigma points
  // Reconstruct mean and covariance matrixes
  PredictMeanAndCovariance();

  PredictRadarMeasurement();

  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

// Perform the appropriate update step depending on sensor type.
void UKF::Update(MeasurementPackage meas_package) {
  if (!use_radar_ && !use_laser_) {
    throw domain_error("UKF::Update: Both use_radar_ and use_laser_ are false. Nothing to update.");
  }

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

  MatrixXd H = MatrixXd(2, 5);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  VectorXd y = z - H * x_;

  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;

  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  x_ = x_ + (K * y);

  MatrixXd I = MatrixXd::Identity(5, 5);

  P_ = (I - K * H) * P_;
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Tc.fill(0);
  VectorXd diff = VectorXd(n_x_);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    diff = Zsig_.col(i) - z_pred_;
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * diff.transpose();
  }

  MatrixXd K = Tc * S_.inverse();
  x_ = x_ + K * (z - z_pred_);
  P_ = P_ - K * S_ * K.transpose();

}



void UKF::GenerateSigmaPoints() {

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig_.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i + 1)     = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig_.col(i + 1+ n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

}


MatrixXd UKF::AugmentSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t) {

  double px, py, v, psi, psi_dot, nu_a, nu_psidd;
  double u_px, u_py, u_v, u_psi, u_psi_dot;
  double n_px, n_py, n_v, n_psi, n_psi_dot;

  for (int i = 0; i < 15; i++) {
    px = Xsig_aug(0, i);
    py = Xsig_aug(1, i);
    v = Xsig_aug(2, i);
    psi = Xsig_aug(3, i);
    psi_dot = Xsig_aug(4, i);
    nu_a = Xsig_aug(5, i);
    nu_psidd = Xsig_aug(6, i);
  
    if (psi_dot != 0) {
      u_px = (v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi));
      u_py = (v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi));
      u_v = 0;
      u_psi = psi_dot * delta_t;
      u_psi_dot = 0;

      n_px = 0.5 * pow(delta_t, 2.0) * cos(psi) * nu_a;
      n_py = 0.5 * pow(delta_t, 2.0) * sin(psi) * nu_a;
      n_v = delta_t * nu_a;
      n_psi = 0.5 * pow(delta_t, 2.0) * nu_psidd;
      n_psi_dot = delta_t * nu_psidd;
    } else {
      u_px = v * cos(psi) * delta_t;
      u_py = v * sin(psi) * delta_t;
      u_v = 0;
      u_psi = psi_dot * delta_t;
      u_psi_dot = 0;

      n_px = 0.5 * pow(delta_t, 2.0) * cos(psi) * nu_a;
      n_py = 0.5 * pow(delta_t, 2.0) * sin(psi) * nu_a;
      n_v = delta_t * nu_a;
      n_psi = 0.5 * pow(delta_t, 2.0) * nu_psidd;
      n_psi_dot = delta_t * nu_psidd;
    }
    
  VectorXd previous = VectorXd(5);
  previous << px, py, v, psi, psi_dot;

  VectorXd update = VectorXd(5);
  update << u_px, u_py, u_v, u_psi, u_psi_dot;

  VectorXd noise = VectorXd(5);
  noise << n_px, n_py, n_v, n_psi, n_psi_dot;

  Xsig_pred_.col(i) = previous + update + noise;

  }
  
}

void UKF::PredictMeanAndCovariance() {
  VectorXd diff = VectorXd(n_x_);
  
  x_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    diff = Xsig_pred_.col(i) - x_;
    P_ += weights_(i) * diff * diff.transpose();
  }  
}


void UKF::PredictRadarMeasurement() {

  double s_px, s_py, s_v, s_psi, s_psi_dot;
  double m_rho, m_psi, m_rho_dot;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    s_px = Xsig_pred_(0, i);
    s_py = Xsig_pred_(1, i);
    s_v = Xsig_pred_(2, i);
    s_psi = Xsig_pred_(3, i);
    s_psi_dot = Xsig_pred_(4, i);

    m_rho = sqrt(pow(s_px, 2.0) + pow(s_py, 2.0));
    m_psi = atan2(s_py, s_px);
    m_rho_dot = ((s_px * cos(s_psi) * s_v) + (s_py * sin(s_psi) * s_v)) / m_rho;

    Zsig_.col(i) << m_rho, m_psi, m_rho_dot;

    z_pred_ += Zsig_.col(i) * weights_(i);
  }

  VectorXd diff = VectorXd(n_z_);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    diff = Zsig_.col(i) - z_pred_;

    while (diff(1) >  M_PI) diff(1) -= 2. * M_PI;
    while (diff(1) < -M_PI) diff(1) += 2. * M_PI;

    S_ += (weights_(i) * diff * diff.transpose());
  }

  S_ += R_;

}

void UKF::UpdateState(VectorXd z) {

}

