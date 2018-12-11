#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() == 0 || ground_truth.size() == 0){
      cout << "ERROR!! Vector size is zero." << endl;
      return rmse;
    }
	if(estimations.size() != ground_truth.size()){
      cout << "ERROR!! Vector size is not equal." << endl;
      return rmse;
    }

    for(int i=0; i < estimations.size(); ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      rmse += residual;
    }

	//Calculate mean and square root of rmse and return it
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}