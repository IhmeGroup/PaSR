#ifndef DROPLETARRAY_H_
#define DROPLETARRAY_H_

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

class DropletArraySolver;

class DropletArray {
  public:
    DropletArray() = default;

    ~DropletArray();

    void initialize(int N, double R0, double k_e, double Re_scale, double c, double omega);

    void computeRHS();

    inline int getSize() const { return this->N; }

    inline double get_N_r2(int i) const { return this->N_r2_arr[i]; }
    inline void set_N_r2(int i, double val) { this->N_r2_arr[i] = val; }
    inline double get_ddt_N_r2(int i) const { return this->ddt_N_r2_arr[i]; }

    void reset_to_N0() { this->N_r2_arr = this->N0_r2_arr; }

    void solve_to_time(double tend);

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

    void calculate_N0();

    void calculate_corr();

    // Inputs
    int N;
    double R0;
    double k_e; // TODO: Dynamically calculate this

    double Re_scale;

    double c;
    double omega;

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