#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_,n_sig_);
  weights_ = VectorXd(n_sig_);
  lambda_ = 3 - n_aug_;
  nis_lidar_ = 0.0;
  nis_radar_ = 0.0;

  P_ << 1, 0, 0, 0, 0,
  		0, 1, 0, 0, 0,
  		0, 0, 1, 0, 0,
  		0, 0, 0, 1, 0,
  		0, 0, 0, 0, 1;

  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
  			  0,std_radphi_ * std_radphi_, 0,
  			  0, 0,std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
  			 0, std_laspy_ * std_laspy_;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

	if(!is_initialized_){


		if (meas_package.sensor_type_== MeasurementPackage::LASER){
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1] , 0, 0, 0;
		}
		else if(meas_package.sensor_type_== MeasurementPackage::RADAR){
			double rho = meas_package.raw_measurements_[0]; 
      		double phi = meas_package.raw_measurements_[1];
      		double rho_dot = meas_package.raw_measurements_[2];

      		double px = rho * cos(phi);
      		double py = rho * sin(phi);

      		//HOW VELOCITY IS CALCULATED HERE
      		
      /*	Note that although radar does include velocity information, 
      		the radar velocity and the CTRV velocity are not the same. 
      		Radar velocity is measured from the autonomous vehicle's perspective. 
      		If you drew a straight line from the vehicle to the bicycle, radar measures the velocity along that line.

			In the CTRV model, the velocity is from the object's perspective, 
			which in this case is the bicycle; the CTRV velocity is tangential to the circle along which the bicycle travels. 
			Therefore, you cannot directly use the radar velocity measurement to initialize the state vector.
      */
      		double vx = rho_dot * cos(phi);
      		double vy = rho_dot * sin(phi);
      		double v = sqrt(vx * vx + vy * vy);

      		x_ << px, py , v, 0, 0;	
		}
		//Saving First timestamp
		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
   *  Prediction
   ****************************************************************************/

	// Calculate the timestamp between measuremetns in seconds

	double dt = (meas_package.timestamp_ - time_us_ )/ 1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(dt);

	/*****************************************************************************
   *  Update
   ****************************************************************************/

	if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
		UpdateRadar(meas_package);
	}
	if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
		UpdateLidar(meas_package);
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // PART 1 : Generation of Sigma Points
	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
	P_aug.fill(0.0);
	MatrixXd Xsig_aug = MatrixXd(n_aug_,n_sig_);
	Xsig_aug.fill(0.0);
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	P_aug.topLeftCorner(5,5) = P_;

	//MatrixXd Q = MatrixXd(2,2);

	//Q << std_a_ * std_a_,0,0,std_yawdd_ * std_yawdd_;

	//P_aug.bottomRightCorner(2,2) = Q;

	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	MatrixXd A = P_aug.llt().matrixL();

	Xsig_aug.col(0) = x_aug;

	for(int i=0;i < n_aug_;++i){
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);

	}

	// PART 2: Prediction of Sigma Points
	//	After Generating the Sigma points, now we have to fill these values into the Process Model to get the predicted Sigma points.
	//	Remember here that Xsig_aug has size 7 X 15 which includes both acc. and yaw noise. When it is passed to Process model, We will get 5 X 15 Xsig_pred_ matrix.
	for(int i=0;i < n_sig_;++i){
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double vel = Xsig_aug(2,i);
		double yaw_rate = Xsig_aug(3,i);
		double change_in_yaw = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_change_in_yaw = Xsig_aug(6,i);
 
		double px_p,py_p;
		if(fabs(change_in_yaw) > 0.001){
			px_p = p_x + vel/change_in_yaw * (sin(yaw_rate + change_in_yaw*delta_t) - sin(yaw_rate));
			py_p = p_y + vel/change_in_yaw * ( cos(yaw_rate) - cos(yaw_rate + change_in_yaw*delta_t) );
		}
		else{
			px_p = p_x + vel * delta_t*cos(yaw_rate);
	        py_p = p_y + vel * delta_t*sin(yaw_rate);
		}

		double vel_p = vel;
		double yaw_p = yaw_rate + change_in_yaw*delta_t;
    	double change_in_yaw_p = change_in_yaw;

    	// Add noise

    	px_p = px_p + 0.5 * nu_a * (delta_t * delta_t) * cos(yaw_rate);
    	py_p = py_p + 0.5 * nu_a * (delta_t * delta_t) * sin(yaw_rate);
    	vel_p = vel_p + delta_t * nu_a;
    	yaw_p = yaw_p + 0.5 * (delta_t * delta_t) * nu_change_in_yaw;
    	change_in_yaw_p = change_in_yaw_p + delta_t * nu_change_in_yaw;

    	Xsig_pred_(0,i) = px_p;
    	Xsig_pred_(1,i) = py_p;
    	Xsig_pred_(2,i) = vel_p;
    	Xsig_pred_(3,i) = yaw_p;
    	Xsig_pred_(4,i) = change_in_yaw_p;

	}

	// PART 3 : Calculate Predicted Mean and Covariance 
	// Calculation of weights

	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < n_sig_; ++i)
	{
		weights_(i) = 0.5 / (lambda_ + n_aug_); 
	}

//	for(int i =0;i<n_sig_;++i){
//      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
//  }

  x_ = Xsig_pred_ * weights_; // 5 X 15 and 15 X 1 --> 5 X 1. Vectorized Summmation

  P_.fill(0.0);
  	for (int i = 0; i < n_sig_; ++i)
  	{
  		VectorXd x_diff = Xsig_pred_.col(i) - x_;

  		//Normalization of angle. As we take difference of state vector, we subtract angle also. Subtracting an angle in Kalman filter
  		//is a problem because result might be 2*pi plus a small angle instead of just a small angle.

  		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    	while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

  		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
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
	MatrixXd zsig = MatrixXd(2,n_sig_);
	//MatrixXd zsig = Xsig_pred_.block(0, 0, 2, n_sig_);
	double px = meas_package.raw_measurements_[0];
	double py = meas_package.raw_measurements_[1];

	z << px,py;

	for (int i = 0; i < n_sig_; ++i)
	{
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);

		zsig(0,i) = px;
		zsig(1,i) = py;
	}

	VectorXd z_pred = VectorXd(2);
	MatrixXd S = MatrixXd(2,2);
	S.fill(0.0);
	z_pred.fill(0.0);
//	for (int i=0; i < n_sig_; i++) {
//      z_pred = z_pred + weights_(i) * zsig.col(i);
//  }
	z_pred = zsig * weights_; // 2 X 15 and 15 X 1 --> 2 X 1. Vectorized Sum

	for(int i=0;i<n_sig_;++i){
      //angle normalization
      VectorXd diff = zsig.col(i) - z_pred;
    	
      S = S + (weights_(i) * diff * diff.transpose());
  }
  S = S + R_lidar_;

  // Calculation of Cross-correlation matrix
  MatrixXd T = MatrixXd(n_x_,2);
  T.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
  	VectorXd x_diff = Xsig_pred_.col(i) - x_;
  	
  	VectorXd z_diff = zsig.col(i) - z_pred;
  	
  	T = T + weights_(i) * x_diff * z_diff.transpose();

  }

  // Calculation of Kalman Gain K
  MatrixXd K_Gain = T * S.inverse();


  //Update State

  x_ = x_ + K_Gain * (z - z_pred);

  //Update Covariance Matrix
  P_ = P_ - K_Gain * S * K_Gain.transpose();

  nis_lidar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

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

	VectorXd z = VectorXd(3);
	double rho = meas_package.raw_measurements_[0];
	double phi = meas_package.raw_measurements_[1];
	double rho_dot = meas_package.raw_measurements_[2];

	z << rho,phi,rho_dot;

	MatrixXd zsig = MatrixXd(3,n_sig_);
	VectorXd z_pred = VectorXd(3);
	MatrixXd S = MatrixXd(3,3);

	for (int i = 0; i < n_sig_; ++i)
	{
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);
		double v = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		zsig(0,i) = sqrt(px * px + py * py);
		zsig(1,i) = atan2(py,px);
		zsig(2,i) = (v * (px * cos(yaw) + py * sin(yaw))) / sqrt(px * px + py * py);


	}
	z_pred.fill(0.0);
	z_pred = zsig * weights_; // 3 X 15 and 15 X 1 --> 3 X 1. Vectorized Sum

	S.fill(0.0);
	for(int i=0;i<n_sig_;++i){
      //angle normalization
      VectorXd diff = zsig.col(i) - z_pred;
    	while (diff(1)> M_PI) {
        	diff(1)-=2.*M_PI;
    	}
    	while (diff(1)<-M_PI) {
    	diff(1)+=2.*M_PI;
		}
      S = S + (weights_(i) * diff * diff.transpose());
  }
  S = S + R_radar_;

  // Calculation of Cross-correlation matrix
  MatrixXd T = MatrixXd(n_x_,3);
  T.fill(0.0);
  for (int i = 0; i < n_sig_; ++i)
  {
  	VectorXd x_diff = Xsig_pred_.col(i) - x_;
  	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

  	VectorXd z_diff = zsig.col(i) - z_pred;
  	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  	
  	T = T + weights_(i) * x_diff * z_diff.transpose();

  }

  // Calculation of Kalman Gain K
  MatrixXd K_Gain = T * S.inverse();


  //Update State

  x_ = x_ + K_Gain * (z - z_pred);

  //Update Covariance Matrix
  P_ = P_ - K_Gain * S * K_Gain.transpose();

  nis_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}
