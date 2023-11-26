#ifndef DROPLETARRAY_H_
#define DROPLETARRAY_H_

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

const int DEFAUlT_N = 1000;
const double DEFAUlT_R0 = 0.0013801623456003516;
const double DEFAUlT_k_e = 0.9/4.;
const double DEFAUlT_Re_scale = 29059.82905982906;
const double DEFAUlT_c = 20.;
const double DEFAUlT_omega = 0.22;
const double DEFAUlT_rho_l = 1000.;

const double DEFAULT_Tsat = 371.53327674824317;
const double DEFAULT_kv = 0.01918528617968421;
const double DEFAULT_rhof = 614.2155646136666;
const double DEFAULT_rhov = 3.4709589577683153;
const double DEFAULT_muv = 7.215094092485614e-06;
const double DEFAULT_hfg = 316884.8797578494;
const double DEFAULT_Tw = 1000.;

class DropletArraySolver;

class DropletArray {
  public:
    DropletArray() = default;

    ~DropletArray();

    void initialize(std::string input_file);

    void initialize(int N, double R0, double k_e, double Re_scale, double c, double omega, double rho_l,
                    double Tsat, double kv, double rhof, double rhov, double muv, double hfg, double Tw);

    void computeRHS();

    inline int getSize() const { return this->N; }

    inline double get_N_r2(int i) const { return this->N_r2_arr[i]; }
    inline void set_N_r2(int i, double val) { this->N_r2_arr[i] = val; }
    inline double get_ddt_N_r2(int i) const { return this->ddt_N_r2_arr[i]; }

    // For testing only
    void solve_to_time(double tend);

    // For interaction with PaSR
    inline double get_initial_mass() {
      return this->rho_l * 4. / 3. * 3.1415 * this->R0 * this->R0 * this->R0;
    }

    void reset_to_N0();
    void take_step(double dt_phys, double Tm);
    int inject(double m_quant);

    void write_output(std::string filename);

  private:
    inline double beta(double r) {
      return std::pow(r, this->c*this->omega) * std::pow(1.0 - r, this->c*(1.0 - this->omega));
    }

    static inline double trapz(Eigen::VectorXd x, Eigen::VectorXd y) {
      double integral = 0.0;
      for (int i = 0; i < x.size() - 1; i++) {
        integral += (x(i+1) - x(i)) * (y(i+1) + y(i)) / 2.0;
      }
      return integral;
    }

    inline double calculate_mass() {
      return this->rho_l * 4. / 3. * 3.1415 * this->R0 * this->R0 * this->R0 
        * DropletArray::trapz(this->r_arr, 
                              this->N_r_arr.array() * this->r3_arr.array());
    }

    void calculate_N0();

    void calculate_corr();

    // Inputs
    int N;
    double R0;
    double k_e; // TODO: Dynamically calculate this

    double Re_scale;

    double c;
    double omega;

    double rho_l;

    double Tsat, kv, rhof, rhov, muv, hfg, Tw;

    double Tm;

    // Compute
    double C;

    Eigen::VectorXd r_arr;
    Eigen::VectorXd r2_arr;
    Eigen::VectorXd r3_arr;

    Eigen::VectorXd N_r_arr;
    Eigen::VectorXd N_r2_arr;

    Eigen::VectorXd corr_1; // Bouyancy
    Eigen::VectorXd corr_2; // Convective
    Eigen::VectorXd ddt_N_r2_arr;


    Eigen::VectorXd N0_r_arr;
    Eigen::VectorXd N0_r2_arr;

    // Solver
    DropletArraySolver* das;

    // Interaction with gas phase
    double m_evap_pool;
    double m_evap_total;
};

class DropletArraySolver {
  public:
    DropletArraySolver() {
      SUNContext_Create(NULL, &this->ctx);
    };

    ~DropletArraySolver() {
      SUNContext_Free(&this->ctx);
      ERKStepFree(&this->erk_mem);
      N_VDestroy_Serial(this->y);
      N_VDestroy_Serial(this->y0);
    };

    void initialize(DropletArray* da) {
      this->da = da;

      // Set up solver
      // Create vector for solution
      this->y = N_VNew_Serial(this->da->getSize(), this->ctx);
      this->y0 = N_VNew_Serial(this->da->getSize(), this->ctx);

      // Set initial conditions
      this->getU(y0);

      // Create solver
      this->erk_mem = ERKStepCreate(&DropletArraySolver::f, 0.0, this->y0, this->ctx);
      ERKStepSetUserData(this->erk_mem, this);

      ERKStepSStolerances(this->erk_mem, 1e-6, 1e-6);
      ERKStepSetOrder(this->erk_mem, 2);
      ERKStepSetMaxNumSteps(this->erk_mem, 100000);
    }

    void solve_to_time(double tend) {
      this->getU(this->y0);

      // Reset solver
      ERKStepReInit(this->erk_mem, &DropletArraySolver::f, 0.0, this->y0);
      ERKStepSetUserData(this->erk_mem, this);

      ERKStepSStolerances(this->erk_mem, 1e-6, 1e-6);
      ERKStepSetOrder(this->erk_mem, 2);
      ERKStepSetMaxNumSteps(this->erk_mem, 100000);

      double t;
      ERKStepEvolve(this->erk_mem, tend, this->y, &t, ARK_NORMAL);

      this->setU(this->y);
      this->da->computeRHS();
    }

  private:
    inline void setU(N_Vector u) {
      realtype*  udata = N_VGetArrayPointer_Serial(u);
      // std::cout << udata[0] << std::endl;
      for (int i = 0; i < this->da->getSize(); ++i) {
        this->da->set_N_r2(i, udata[i]);
      }
      this->da->computeRHS();
    }

    inline void getU(N_Vector u) {
      realtype*  udata = N_VGetArrayPointer_Serial(u);
      for (int i = 0; i < this->da->getSize(); ++i) {
        udata[i] = this->da->get_N_r2(i);
      }
    }

    inline void getdU(N_Vector du) {
      realtype*  udata = N_VGetArrayPointer_Serial(du);
      for (int i = 0; i < this->da->getSize(); ++i) {
        udata[i] = this->da->get_ddt_N_r2(i);
      }
    }

    static int f(realtype t, N_Vector u, N_Vector du, void* user_data) {
      DropletArraySolver* das = static_cast<DropletArraySolver*>(user_data);

      das->setU(u);
      das->getdU(du);

      return 0;
    }

    DropletArray* da;

    // ARKode
    SUNContext ctx;
    void* erk_mem;

    N_Vector y0;
    N_Vector y;
};

#endif /* DROPLETARRAY_H_ */