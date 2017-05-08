#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  MatrixXd Xsig_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  MatrixXd Zsig_;

  VectorXd z_pred_;

  ///* state covariance matrix
  MatrixXd P_;

  MatrixXd R_;

  MatrixXd S_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  //Measurement dimension
  int n_z_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  long long previous_timestamp_;

  bool initialized_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);


  /**
   * Extracts measurement data depending on sensor type, updates delta_t, and returns
   * the measurement vector.
   * @param meas_package The measurement at k+1
   * @param delta_t where the calculated elapsed time will be stored
   */
  VectorXd ExtractMeasurementData(MeasurementPackage meas_package, long long &previous_timestamp, double &delta_t);

  /**
   * Samples sigma points to be fed through the process model
   */
  void GenerateSigmaPoints();

  /**
   * Augments the sigma points with process error
   */
  MatrixXd AugmentSigmaPoints();

  void PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t);

  void PredictMeanAndCovariance();

  void PredictRadarMeasurement();
  
  /**
   * Feed sigma points through the process model
   */
  // void SigmaPointPrediction(MatrixXd* Xsig_out);

  /**
   * Calls UpdateLidar or UpdateRadar depending on sensor type
   * @param meas_package The measurement at k+1
   */
  void Update(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void UpdateState(VectorXd z);
};

#endif /* UKF_H */