#include <iostream>
#include <string>
#include <vector>

#include "cantera/base/Solution.h"

enum MixingModel {NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, EMST};

class Solver {
public:
    explicit Solver(const std::string& input_filename_);
    ~Solver();
    void run();
    void print();

protected:
    void parseInput();
    void takeStep();
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
    unsigned int n_particles;
    int n_steps;
    double t_stop;
    double dt;
    double t;
    unsigned int step;
    MixingModel mixing_model;
    double P;
    std::string comp_fuel, comp_ox;
    double T_fuel, T_ox;
    double phi_global;
    double tau_res, tau_mix;

    Cantera::Solution* gas = nullptr;
    std::vector<double> data;

    unsigned int n_threads;

private:
    
};