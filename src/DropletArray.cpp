#include "DropletArray.h"

#include <fstream>
#include <iomanip>

DropletArray::~DropletArray() {
  delete this->das;
}

//--------------------------------------------------------------------------

void DropletArray::initialize(int N, double R0, double k_e, double Re_scale, double c, double omega) {
  this->N = N;
  this->R0 = R0;
  this->k_e = k_e;
  this->Re_scale = Re_scale;

  this->c = c;
  this->omega = omega;

  this->r2_arr = Eigen::VectorXd::LinSpaced(this->N, 1e-12, 1.0);
  this->r_arr = this->r2_arr.array().sqrt();
  this->r3_arr = this->r_arr.array().cube();

  this->N_r_arr = Eigen::VectorXd::Zero(this->N);
  this->N_r2_arr = Eigen::VectorXd::Zero(this->N);

  this->corr_1 = Eigen::VectorXd::Zero(this->N);
  this->corr_2 = Eigen::VectorXd::Zero(this->N);
  this->ddt_N_r2_arr = Eigen::VectorXd::Zero(this->N);

  this->N0_r_arr = Eigen::VectorXd::Zero(this->N);
  this->N0_r2_arr = Eigen::VectorXd::Zero(this->N);

  this->calculate_N0();
  this->calculate_corr();

  this->das = new DropletArraySolver();
  this->das->initialize(this);
}

//--------------------------------------------------------------------------

void DropletArray::calculate_N0() {
  for (int i = 0; i < this->N; ++i) {
    this->N0_r_arr[i] = this->beta(this->r_arr[i]);
  }

  this->C = 1. / DropletArray::trapz(this->r_arr, this->N0_r_arr.array() * this->r3_arr.array());

  this->N0_r_arr *= this->C;

  this->N0_r2_arr = this->N0_r_arr.array() / 2. / this->r_arr.array();

}

//--------------------------------------------------------------------------

void DropletArray::calculate_corr() {
  for (int i = 0; i < this->N; ++i) {
    double r_actual = this->r_arr[i] * this->R0;

    this->corr_1[i] = 1. + 0.276 * std::pow(this->Re_scale * r_actual, 0.5);
    this->corr_2[i] = 1. + 8.5 * std::pow(2. * r_actual * 100., 3.*0.45);
  }
}

//--------------------------------------------------------------------------

void DropletArray::computeRHS() {
  // Calculate N
  N_r_arr = N_r2_arr.array() * 2. * this->r_arr.array();

  // Calculate dN/dt
  for (int i = 0; i < N; ++i) {
    double FluxL, FluxR;
    if (i < this->N-1)
      FluxR = -this->N_r2_arr[i+1] * this->corr_1[i+1] * this->corr_2[i+1];
    else 
      FluxR = 0.;

    FluxL = -this->N_r2_arr[i] * this->corr_1[i] * this->corr_2[i];

    double h;
    if (i > 0 && i < N-1)
      h = (this->r2_arr[i+1] - this->r2_arr[i-1]) / 2.;
    else if (i == 0)
      h = this->r2_arr[i+1] - this->r2_arr[i];
    else if (i == N-1)
      h = this->r2_arr[i] - this->r2_arr[i-1];

    ddt_N_r2_arr[i] = -(FluxR - FluxL) / h;
  }
}

//--------------------------------------------------------------------------

void DropletArray::solve_to_time(double tend) {
  this->reset_to_N0();
  this->computeRHS();

  this->das->solve_to_time(tend);

  std::cout << N_r2_arr << std::endl;

  this->write_output("test.csv");
}

//--------------------------------------------------------------------------

void DropletArray::write_output(std::string filename) {
  std::ofstream myfile;
  myfile.open(filename);

  // Preamble
  myfile << "r,r2,N_r,N_r2" << std::endl;

  for (int i = 0; i < this->N; ++i) {
    myfile << std::setprecision(9) << this->r_arr[i] << ","
            << std::setprecision(9) << this->r2_arr[i] << ","
            << std::setprecision(9) << this->N_r_arr[i] << ","
            << std::setprecision(9) << this->N_r2_arr[i] << std::endl;
  }

  myfile.close();
}






