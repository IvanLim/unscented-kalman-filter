#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // Sanity checks
  int size_est = estimations.size();
  int size_gnd = ground_truth.size();

  int vector_size = estimations[0].size();
	VectorXd rmse(vector_size);
	rmse.fill(0);

  if (size_est == 0 || size_gnd == 0 || size_est != size_gnd) {
      cout << "Tools::CalculateRMSE Error: Estimation or Ground Truth vectors invalid." << endl;
      return rmse;
  }

  //accumulate squared residuals
  VectorXd diff(vector_size);
  for(int i = 0; i < size_est; ++i) {
      diff = estimations[i] - ground_truth[i];
        rmse = rmse.array() + (diff.array() * diff.array());
  }

  //calculate the mean and square root
  rmse = rmse.array() / size_est;
  rmse = rmse.array().sqrt();	

  return rmse;
}
