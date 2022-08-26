#include <iostream>
#include <string>
#include <vector>

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

protected:
    void parseInput();
    void takeStep();
    void subStepInflow();
    void subStepMix();
    void subStepReact();
    void incrementAge();
    void recycleParticle(const unsigned int& ip, const double& p_inj);
    std::vector<double> favreMeanState();
    std::vector<double> meanState();
    bool runDone();

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
    double t_stop;
    double dt;
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

    std::vector<unsigned int>seedvec;
    double p_out, p_mix;

private:
    
};
