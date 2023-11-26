#include "DropletArray.h"

#include <fstream>
#include <iomanip>

#include "../cpptoml/include/cpptoml.h"

DropletArray::~DropletArray() {
  delete this->das;
}

//--------------------------------------------------------------------------

void DropletArray::initialize(std::string input_file) {
  auto config = cpptoml::parse_file(input_file);
  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << "DropletArray parameters:" << std::endl;
  int N = config->get_qualified_as<int>("DropletArray.N").value_or(DEFAUlT_N);
  std::cout << "> DropletArray.N = " << N << std::endl;
  double R0 = config->get_qualified_as<double>("DropletArray.R0").value_or(DEFAUlT_R0);
  std::cout << "> DropletArray.R0 = " << R0 << std::endl;
  double k_e = config->get_qualified_as<double>("DropletArray.k_e").value_or(DEFAUlT_k_e);
  std::cout << "> DropletArray.k_e = " << k_e << std::endl;
  double Re_scale = config->get_qualified_as<double>("DropletArray.Re_scale").value_or(DEFAUlT_Re_scale);
  std::cout << "> DropletArray.Re_scale = " << Re_scale << std::endl;
  double c = config->get_qualified_as<double>("DropletArray.c").value_or(DEFAUlT_c);
  std::cout << "> DropletArray.c = " << c << std::endl;
  double omega = config->get_qualified_as<double>("DropletArray.omega").value_or(DEFAUlT_omega);
  std::cout << "> DropletArray.omega = " << omega << std::endl;
  double rho_l = config->get_qualified_as<double>("DropletArray.rho_l").value_or(DEFAUlT_rho_l);
  std::cout << "> DropletArray.rho_l = " << rho_l << std::endl;
  
  double  Tsat = config->get_qualified_as<double>("DropletArray.Tsat").value_or(DEFAULT_Tsat);
  std::cout << "> DropletArray.Tsat = " << Tsat << std::endl;
  double kv = config->get_qualified_as<double>("DropletArray.kv").value_or(DEFAULT_kv);
  std::cout << "> DropletArray.kv = " << kv << std::endl;
  double rhof = config->get_qualified_as<double>("DropletArray.rhof").value_or(DEFAULT_rhof);
  std::cout << "> DropletArray.rhof = " << rhof << std::endl;
  double rhov = config->get_qualified_as<double>("DropletArray.rhov").value_or(DEFAULT_rhov);
  std::cout << "> DropletArray.rhov = " << rhov << std::endl;
  double muv = config->get_qualified_as<double>("DropletArray.muv").value_or(DEFAULT_muv);
  std::cout <<  "> DropletArray.muv = " << muv << std::endl;
  double hfg = config->get_qualified_as<double>("DropletArray.hfg").value_or(DEFAULT_hfg);
  std::cout << "> DropletArray.hfg = " << hfg << std::endl;
  double Tw = config->get_qualified_as<double>("DropletArray.Ts").value_or(DEFAULT_Tw);
  std::cout << "> DropletArray.Ts = " << Tw << std::endl;

  this->initialize(N, R0, k_e, Re_scale, c, omega, rho_l,
                   Tsat, kv, rhof, rhov, muv, hfg, Tw);
}

//--------------------------------------------------------------------------

void DropletArray::initialize(int N, double R0, double k_e, double Re_scale, double c, double omega, double rho_l,
                              double Tsat, double kv, double rhof, double rhov, double muv, double hfg, double Tw) {
  this->N = N;
  this->R0 = R0;
  this->k_e = k_e;
  this->Re_scale = Re_scale;

  this->c = c;
  this->omega = omega;

  this->rho_l = rho_l;

  this->Tsat = Tsat;
  this->kv = kv;
  this->rhof = rhof;
  this->rhov = rhov;
  this->muv = muv;
  this->hfg = hfg;
  this->Tw = Tw;

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

  this->reset_to_N0();

  this->m_evap_pool = 0.;
  this->m_evap_total = 0.;

  this->Tm = 300.;
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
  double A1 = kv * (Tw - Tsat) / rhof / hfg;
  double A2 = rhov * 9.81 / muv;
  for (int i = 0; i < this->N; ++i) {
    double r_actual = this->r_arr[i] * this->R0;
    double d_drop = 2. * r_actual;

    double Dmax_D0 = 1.;
    double Dmax = Dmax_D0 * d_drop;

    double C = 0.0024;
    double f = 2. * C * d_drop / A1 + 3.76 / std::pow(A1,0.75) / std::pow(A2,0.25) * 5. / 9. * std::pow( std::pow(d_drop/2.,2.25) - 9./16.*C*Dmax/2.*std::pow(A2,0.25)/std::pow(A1,0.25)*std::pow(d_drop,2),-4./9.) * (2.25/std::pow(2.,2.25)*std::pow(d_drop,1.25) - 9./16.*C*std::pow(A2/A1,0.25)*Dmax_D0*3./2.*std::pow(d_drop,2));

    this->corr_1[i] = r_actual / f / (this->k_e*1e-6);



    this->corr_2[i] = (1. + 0.276 * std::pow(this->Re_scale * r_actual, 0.5)) * (1. + 8.5 * std::pow(2. * r_actual * 100., 3.*0.45)) / 2.;
  }
}

//--------------------------------------------------------------------------

void DropletArray::computeRHS() {
  // Calculate N
  N_r_arr = N_r2_arr.array() * 2. * this->r_arr.array();

  double scaleEvap = std::pow(this->Tm / 2200., 3.);

  // Calculate dN/dt
  for (int i = 0; i < N; ++i) {
    double FluxL, FluxR;
    if (i < this->N-1)
      FluxR = -this->N_r2_arr[i+1] * (this->corr_1[i+1] + scaleEvap*this->corr_2[i+1]);
    else 
      FluxR = 0.;

    FluxL = -this->N_r2_arr[i] * (this->corr_1[i] + scaleEvap*this->corr_2[i]);

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

  this->das->solve_to_time(tend);

  std::cout << this->corr_1 << std::endl;

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

//--------------------------------------------------------------------------

void DropletArray::reset_to_N0() {
  this->N_r2_arr = this->N0_r2_arr;
  this->computeRHS();
}

//--------------------------------------------------------------------------

void DropletArray::take_step(double dt_phys, double Tm) {
  this->Tm = Tm;
  // std::cout << "Tm = " << Tm << std::endl;
  // Convert physical time step to dimensionless time step
  double dt = dt_phys / (this->R0 * this->R0 / this->k_e * 1e6);
  double m_pre = this->calculate_mass();
  this->das->solve_to_time(dt);
  double m_post = this->calculate_mass();

  // Add to m_evap
  this->m_evap_pool -= m_post - m_pre;
}

int DropletArray::inject(double m_quant) {
  int nInject = (int) (this->m_evap_pool / m_quant);
  this->m_evap_pool -= nInject * m_quant;

  this->m_evap_total += nInject * m_quant;

  return nInject;
}





