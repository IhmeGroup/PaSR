#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include "cantera/base/ctexceptions.h"
#include "cantera/core.h"

#include "Injector.h"
#include "Particle.h"

enum MixingModel {NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, EMST};

class PartiallyStirredReactor {
public:
    explicit PartiallyStirredReactor(const std::string& input_filename_);
    ~PartiallyStirredReactor();
    void initialize();
    void run();
    void print();
    void check();

protected:
    void parseInput();
    void takeStep();
    void subStepInflow(double dt);
    void subStepMix(double dt);
    void subStepReact(double dt);
    void incrementAge();
    void recycleParticle(unsigned int ip, double p_inj);
    bool runDone();

    double mean(std::function<double(int)> xfunc, bool favre=true);
    double meanState(int iv, bool favre=true);
    void meanState(std::vector<double>* xsumvec, bool favre=true);

    double min(std::function<double(int)> xfunc);
    double minState(int iv);
    void minState(std::vector<double>* minvec);

    double max(std::function<double(int)> xfunc);
    double maxState(int iv);
    void maxState(std::vector<double>* maxvec);

    double meanAge(bool favre=true);
    double minAge();
    double maxAge();

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

    std::vector<unsigned int>seedvec;
    double p_out, p_mix;

    unsigned int check_interval;

private:
    
};
