#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.0;

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
  
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  
  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  // Log files
  logfile_laser_.open("laser_nis.txt", std::ios_base::app);
  logfile_radar_.open("radar_nis.txt", std::ios_base::app);
    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
  
    x_ = VectorXd(5);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      x_(0) = rho * cos(theta);
      x_(1) = rho * sin(theta);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    
	  P_ = MatrixXd(5, 5);
	  P_ << 1, 0, 0, 0, 0,
			    0, 1, 0, 0, 0,
			    0, 0, 1, 0, 0,
			    0, 0, 0, 1, 0,
			    0, 0, 0, 0, 1;
		
    time_us_ = meas_package.timestamp_;
    
    cout << "Initialized, x: " << x_ << endl;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/   
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_) {
      UpdateRadar(meas_package);
    }
  }
  else {
    if (use_laser_) {
      UpdateLidar(meas_package);
    }
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl << endl;
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
  // Generate sigma points
  MatrixXd Xsig_aug;
  AugmentedSigmaPoints(&Xsig_aug);
  
  // Predict sigma points
  SigmaPointPrediction(Xsig_aug, delta_t);
  
  // Predict state x_ and state covariance matrix P_
  PredictMeanAndCovariance();//&x_, &P_);
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;
  unsigned int n_z = 2;
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_, 0,
       0,                     std_laspy_*std_laspy_;
  
  PredictSensorMeasurement(&Zsig, &z_pred, &S, n_z, 1, R);
  UpdateState(Zsig, z_pred, S, meas_package.raw_measurements_, n_z);
  
  double nis = CalculateNIS(z_pred, S, meas_package.raw_measurements_);
  //laser_nis_.push_back(nis);
  logfile_laser_ << nis << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;
  unsigned int n_z = 3;
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0,                       0,
       0,                   std_radphi_*std_radphi_, 0,
       0,                   0,                       std_radrd_*std_radrd_;  
       
  PredictSensorMeasurement(&Zsig, &z_pred, &S, n_z, 0, R);
  UpdateState(Zsig, z_pred, S, meas_package.raw_measurements_, n_z);
  
  double nis = CalculateNIS(z_pred, S, meas_package.raw_measurements_);
  //radar_nis_.push_back(nis);
  logfile_radar_ << nis << endl;
}


/**
 * Prediction functions
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  double sqrt_val = sqrt(lambda_+n_aug_);
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt_val * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_val * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;

}


void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug, double delta_t) {
  
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}


void UKF::PredictMeanAndCovariance() {
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  // It might block here if the x_diff(3) angle is too large, as the angle wrap function is pretty stupid...
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    
    //angle normalization    
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }
  
  //write result
  x_ = x;
  P_ = P;
}


/**
 * Laser and radar update function
 */

void UKF::PredictSensorMeasurement(MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out, unsigned int n_z, unsigned int sensor, MatrixXd& R) {
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    // measurement model
    if (sensor == 0) {                           // Radar
      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;
      double rho = sqrt(p_x*p_x + p_y*p_y);
      Zsig(0,i) = rho;                           //rho
      Zsig(1,i) = atan2(p_y,p_x);                //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2) / rho;       //rho_dot
    }
    else {                                       // Laser
      Zsig(0,i) = p_x;                           // px
      Zsig(1,i) = p_y;                           // py
    }
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R;
  
  //write result
  *Zsig_out = Zsig;
  *z_out = z_pred;
  *S_out = S;
}


void UKF::UpdateState(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S, VectorXd& z, unsigned int n_z) {
    
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
}


/**
 * Calculate NIS
 */
double UKF::CalculateNIS(VectorXd& z_pred, MatrixXd& S, VectorXd& z) {
  VectorXd z_diff = z - z_pred;
  return (z_diff.transpose() * S.inverse() * z_diff);

}

