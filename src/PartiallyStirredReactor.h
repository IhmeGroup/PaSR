#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include "cantera/base/ctexceptions.h"
#include "cantera/core.h"

#include "Injector.h"
#include "Particle.h"

const std::string DEFAULT_MECH_NAME = "";
const unsigned int DEFAULT_N_PARTICLES = 100;
const unsigned int DEFAULT_N_STEPS = -1;
const double DEFAULT_T_STOP = -1.0;
const double DEFAULT_DT = -1.0;
const double DEFAULT_DT_SUB = -1.0;
const std::string DEFAULT_CONVERGENCE_METRIC = "MEAN";
const double DEFAULT_RTOL = -1.0;
const std::string DEFAULT_MIX_MODEL = "FULL_MIX";
const double DEFAULT_PRESSURE = 101325.0;
const std::string DEFAULT_COMP_FUEL = "";
const std::string DEFAULT_COMP_OX = "";
const double DEFAULT_T_FUEL = 300.0;
const double DEFAULT_T_OX = 300.0;
const double DEFAULT_T_INIT = 300.0;
const double DEFAULT_PHI_GLOBAL = 1.0;
const double DEFAULT_TAU_RES = 1.0;
const double DEFAULT_TAU_MIX = 1.0;
const unsigned int DEFAULT_CHECK_INTERVAL = 1;
const bool DEFAULT_CHECK_VERBOSE = false;

enum MixingModel {NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, EMST};
enum ConvergenceMetric {MEAN, MEAN_VAR, HIST};

class PartiallyStirredReactor {
public:
    explicit PartiallyStirredReactor(const std::string& input_filename_);
    ~PartiallyStirredReactor();
    void initialize();
    void run();
    void print();
    void check();

    std::string varName(int iv);

    double minState(int iv);
    void minState(std::vector<double>* minvec);
    double maxState(int iv);
    void maxState(std::vector<double>* maxvec);
    double meanState(int iv, bool favre=true);
    void meanState(std::vector<double>* xsumvec, bool favre=true);
    double varState(int iv, bool favre=true);
    void varState(std::vector<double>* xsumvec, bool favre=true);
    void histState(std::vector<double>* histvec, int iv);
    void histState(std::vector<std::vector<double>>* histvec);

    double meanAge(bool favre=true);
    double minAge();
    double maxAge();

    double meanT(bool favre=true);
    double minT();
    double maxT();

    double meanZ(bool favre=true);
    double minZ();
    double maxZ();

protected:
    void parseInput();
    void takeStep();
    void subStepInflow(double dt);
    void subStepMix(double dt);
    void subStepReact(double dt);
    void incrementAge();
    void recycleParticle(unsigned int ip, double p_inj);
    void calcConvergence();
    bool runDone();

    double min(std::function<double(int)> xfunc);
    double min(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc);

    double max(std::function<double(int)> xfunc);
    double max(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc);

    double mean(std::function<double(int)> xfunc, bool favre=true);
    double mean(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool favre);

    double var(std::function<double(int)> xfunc, bool favre=true);
    double var(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool favre);

    void hist(std::vector<double>* histvec, std::function<double(int)> xfunc);
    void hist(std::vector<double>* histvec, std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc);

    std::string convergenceMetricString(ConvergenceMetric convergence_metric_) {
        switch(convergence_metric_) {
            case MEAN: return "MEAN";
            case MEAN_VAR: return "MEAN_VAR";
            case HIST: return "HIST";
        }
    }

    std::string mixingModelString(MixingModel mixing_model_) {
        switch(mixing_model_) {
            case NO_MIX: return "NO_MIX";
            case FULL_MIX: return "FULL_MIX";
            case CURL: return "CURL";
            case MOD_CURL: return "MOD_CURL";
            case IEM: return "IEM";
            case EMST: return "EMST";
        }
    }

    std::string input_filename;
    std::string mech_filename;
    unsigned int np;
    int n_steps;
    int n_substeps;
    double t_stop;
    double dt_step;
    double dt_sub_target;
    ConvergenceMetric convergence_metric;
    double rtol, rerror;
    double t;
    unsigned int step;
    MixingModel mixing_model;
    double P;
    std::string comp_fuel, comp_ox;
    double T_fuel, T_ox, T_init;
    double h_fuel, h_ox;
    double phi_global;
    double tau_res, tau_mix;

    std::vector<std::shared_ptr<Cantera::Solution>> solvec;
    std::vector<std::shared_ptr<Cantera::ThermoPhase>> gasvec;
    std::vector<Cantera::IdealGasConstPressureReactor*> reactorvec;
    std::vector<Cantera::ReactorNet*> rnetvec;

    std::vector<double> xvec;
    std::vector<Particle> pvec;
    std::vector<Injector> injvec;
    unsigned int nsp;
    unsigned int nv;
    std::vector<double> Y_fuel, Y_ox, Y_phi;
    std::vector<double> xtemp;
    std::vector<double> xmean_old;

    std::vector<unsigned int> i_fuel;
    std::vector<unsigned int> iv_check;
    std::vector<unsigned int> seedvec;
    double p_out, p_mix;

    unsigned int check_interval;
    bool check_verbose;

private:
    
};
