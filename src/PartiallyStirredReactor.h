#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <random>
#include <unordered_map>

#include "cantera/base/ctexceptions.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ReactorNet.h"

#include "../cpptoml/include/cpptoml.h"

#include "Histogram.h"
#include "Injector.h"
#include "Particle.h"

const std::string DEFAULT_MECH_NAME = "";
const bool DEFAULT_RESTART = false;
const std::string DEFAULT_RESTART_NAME = "";
const double DEFAULT_P_STOICH = 0.6;
const int DEFAULT_N_PARTICLES = 100;
const int DEFAULT_N_STEPS = -1;
const double DEFAULT_T_STOP = -1.0;
const double DEFAULT_T_EXTINCT = -1.0;
const std::string DEFAULT_T_EXTINCT_MODE = "MEAN";
const double DEFAULT_DT = -1.0;
const double DEFAULT_DT_SUB = -1.0;
const std::string DEFAULT_CONVERGENCE_METRIC = "MEAN";
const int DEFAULT_STATS_WINDOW = 1;
const double DEFAULT_RTOL = -1.0;
const int DEFAULT_MIN_STEPS_CONVERGE = -1;
const std::string DEFAULT_MIX_MODEL = "FULL_MIX";
const std::string DEFAULT_INJECTION_MODE = "NONPREMIXED";
const double DEFAULT_PILOT_FLOW = 0.0;
const double DEFAULT_PILOT_T_STOP = -1.0;
const double DEFAULT_PRESSURE = 101325.0;
const std::string DEFAULT_COMP_FUEL = "";
const std::string DEFAULT_COMP_OX = "";
const double DEFAULT_T_FUEL = 300.0;
const double DEFAULT_T_OX = 300.0;
const double DEFAULT_T_INIT = 300.0;
const double DEFAULT_PHI_GLOBAL = 1.0;
const std::string DEFAULT_TAU_RES_MODE = "EXP_MEAN";
const double DEFAULT_TAU_RES_VALUE = 1.0;
const std::string DEFAULT_TAU_RES_HIST_NAME = "";
const double DEFAULT_TAU_MIX = 1.0;
const unsigned int DEFAULT_CHECK_INTERVAL = 1;
const std::vector<std::string> DEFAULT_CHECK_VARIABLE_NAMES{};
const bool DEFAULT_CHECK_VERBOSE = false;
const int DEFAULT_WRITE_INTERVAL = -1;

enum class MixingModel {NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, EMST_1D, EMST, KER_M};
enum class InjectionMode {PREMIXED, NONPREMIXED};
enum class ConvergenceMetric {MEAN, MEAN_VAR, HIST};
enum class ExtinctMode {MEAN, MAX, _COUNT};
enum class TauResMode {EXP_MEAN, DISTRIBUTION};
enum class VariableType {STATE, AUX, DERIVED};

const std::unordered_map<std::string, MixingModel> mixing_model_map = {
    {"NO_MIX", MixingModel::NO_MIX},
    {"FULL_MIX", MixingModel::FULL_MIX},
    {"CURL", MixingModel::CURL},
    {"MOD_CURL", MixingModel::MOD_CURL},
    {"IEM", MixingModel::IEM},
    {"EMST_1D", MixingModel::EMST_1D},
    {"EMST", MixingModel::EMST},
    {"KER_M", MixingModel::KER_M}
};

const std::unordered_map<std::string, InjectionMode> injection_mode_map = {
    {"PREMIXED", InjectionMode::PREMIXED},
    {"NONPREMIXED", InjectionMode::NONPREMIXED}
};

const std::unordered_map<std::string, ConvergenceMetric> convergence_metric_map = {
    {"MEAN", ConvergenceMetric::MEAN},
    {"MEAN_VAR", ConvergenceMetric::MEAN_VAR},
    {"HIST", ConvergenceMetric::HIST}
};

const std::unordered_map<std::string, ExtinctMode> extinct_mode_map = {
    {"MEAN", ExtinctMode::MEAN},
    {"MAX", ExtinctMode::MAX}
};

const std::unordered_map<std::string, TauResMode> tau_res_mode_map = {
    {"EXP_MEAN", TauResMode::EXP_MEAN},
    {"DISTRIBUTION", TauResMode::DISTRIBUTION}
};

const std::unordered_map<std::string, VariableType> variable_type_map = {
    {"STATE", VariableType::STATE},
    {"AUX", VariableType::AUX},
    {"DERIVED", VariableType::DERIVED}
};

const std::string RAW_NAME = "particle_data";
const std::string RAW_EXT = ".csv";
const std::string STATS_DIR = "stats";
const std::string STATS_PREFIX = "stats_";
const std::string STATS_EXT = ".csv";

const int WRITE_PRECISION = 15;

#pragma omp declare \
    reduction( \
        vec_particle_plus : \
        std::vector<Particle> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            std::plus<Particle>())) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))

class PartiallyStirredReactor {
public:
    explicit PartiallyStirredReactor(const std::string& input_filename_);
    ~PartiallyStirredReactor();
    void initialize();
    void run();
    void print();
    void check(bool force=false);

    std::string variableName(int iv);
    int variableIndex(std::string name);
    int nVariables() { return n_state_variables + n_aux_variables + n_derived_variables; }

    double min(int iv, bool all=false);
    double max(int iv, bool all=false);
    double mean(int iv, bool all=false, bool favre=true);
    double variance(int iv, bool all=false, bool favre=true);
    double variance(int iv, double meanval, bool all=false);
    void hist(std::vector<double>* histvec, int iv, bool all=false);

    void minState(std::vector<double>* minvec, bool all=false);
    void maxState(std::vector<double>* maxvec, bool all=false);
    void meanState(std::vector<double>* xsumvec, bool all=false, bool favre=true);
    void meanState(std::vector<Particle>* pvec_, std::vector<double>* xsumvec, bool favre=true);
    void varianceState(std::vector<double>* xvarvec, bool all=false, bool favre=true);
    void varianceState(std::vector<Particle>* pvec_, std::vector<double>* xvarvec, bool favre=true);
    void varianceState(std::vector<double>* xvarvec, std::vector<double>* xmeanvec, bool all=false);
    void histState(std::vector<std::vector<double>>* histvec, bool all=false);

protected:
    void parseInput();
    template <typename T> T parseValue(
        std::shared_ptr<cpptoml::table> input,
        std::string name,
        T default_value);
    template <typename T> T parseEnumValue(
        std::shared_ptr<cpptoml::table> input,
        std::string name,
        std::unordered_map<std::string, T> enum_map,
        std::string default_value_str);
    template <typename T> std::vector<T> parseArray(
        std::shared_ptr<cpptoml::table> input,
        std::string name,
        std::vector<T> default_value);
    void readRestart();
    void takeStep();
    void calcDt();
    void subStepInflow(double dt);
    void subStepMix(double dt);
    void subStepReact(double dt);
    void incrementAge();
    void recycleParticle(unsigned int ip, double p_inj, int tid=0);
    void calcConvergence();
    bool runDone();
    void copyState();
    std::string statsPath(std::string name);
    void writeStatsHeaders();
    void writeStats(bool force=false);
    void writeRawHeaders();
    void writeRaw(bool force=false);

    void addVariable(
        std::string name,
        std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> getter,
        VariableType type);
    void addAuxVariable(
        std::string name,
        std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> getter);
    void addDerivedVariable(
        std::string name,
        std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> getter);
    void addCheckVariable(int iv);
    void checkVariable(int iv);

    double min(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all=false);
    double max(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all=false);
    double mean(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all=false, bool favre=true);
    double variance(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all=false, bool favre=true);
    void hist(std::vector<double>* histvec, std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all=false);

    std::string input_filename;
    std::string mech_filename;
    bool restart;
    std::string restart_filename;
    double p_phi_equil;
    unsigned int n_particles;
    int n_steps;
    int n_substeps;
    int min_steps_converge;
    double t_stop;
    double T_extinct;
    ExtinctMode T_extinct_mode;
    double dt_step;
    int n_sub;
    double dt_sub, dt_sub_target;
    ConvergenceMetric convergence_metric;
    int n_stat;
    int i_stat;
    double rtol, rerror;
    double t;
    unsigned int step;
    MixingModel mixing_model;
    InjectionMode injection_mode;
    double pilot_flow;
    double pilot_t_stop;
    double P;
    std::string comp_fuel, comp_ox;
    double T_fuel, T_ox, T_mix;
    double h_fuel, h_ox, h_mix;
    double phi_global;
    TauResMode tau_res_mode;
    double tau_res_value;
    std::string tau_res_hist_name;
    Histogram tau_res_hist;
    double tau_mix;

    bool run_done;
    bool particles_injected;
    unsigned int n_recycled, n_recycled_check;

    std::vector<std::shared_ptr<Cantera::Solution>> solvec;
    std::vector<std::shared_ptr<Cantera::ThermoPhase>> gasvec;
    std::vector<Cantera::IdealGasConstPressureReactor*> reactorvec;
    std::vector<Cantera::ReactorNet*> rnetvec;

    std::vector<Particle> pvec, pvec_temp1, pvec_temp2, pvec_partemp;
    std::vector<Injector> injvec;
    unsigned int n_species;
    unsigned int n_state_variables, n_aux_variables, n_derived_variables;
    unsigned int id_iterator;
    std::vector<unsigned int> i_fuel;
    std::vector<double> Y_fuel, Y_ox, Y_mix;
    std::vector<double> xtemp1, xtemp2;
    std::vector<double> xmean_old, xvar_old;

    std::vector<std::string> aux_variable_names, derived_variable_names;
    std::vector<std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)>> variable_functions;
    std::vector<double> xtemp_derived;

    std::vector<std::uniform_int_distribution<unsigned int>> dists_uni_int;
    std::vector<std::uniform_real_distribution<double>> dists_uni_real;
    std::vector<std::mt19937> rand_engines;
    double p_out, p_mix;

    unsigned int check_interval;
    std::vector<std::string> check_variable_names;
    std::vector<unsigned int> iv_check;
    bool check_verbose;
    unsigned int write_raw_interval;
    unsigned int write_stats_interval;

private:
    
};