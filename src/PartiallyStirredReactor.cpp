#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include "../cpptoml/include/cpptoml.h"

#include "common.h"
#include "PartiallyStirredReactor.h"

#ifdef __APPLE__
#include <filesystem>
namespace fs = std::__fs::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

PartiallyStirredReactor::PartiallyStirredReactor(const std::string& input_filename_) :
    input_filename(input_filename_), run_done(false),
    step(0), t(0.0), p_out(0.0), i_stat(1), id_iterator(0), n_particles(0),
    n_state_variables(0), n_aux_variables(0), n_derived_variables(0), n_species(0),
    particles_injected(false), n_recycled(0), n_recycled_check(0)
{
    parseInput();
}

void PartiallyStirredReactor::parseInput() {
    auto config = cpptoml::parse_file(input_filename);

    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Input file parameters:" << std::endl;

    // Mechanism
    mech_filename = config->get_qualified_as<std::string>("Mechanism.name").value_or(DEFAULT_MECH_NAME);
    std::cout << "> Mechanism.name = " << mech_filename << std::endl;

    // Initialization
    restart = config->get_qualified_as<bool>("Initialization.restart").value_or(DEFAULT_RESTART);
    std::cout << "> Initialization.restart = " << restart << std::endl;
    if (restart) {
        restart_filename = config->get_qualified_as<std::string>("Initialization.name").value_or(DEFAULT_RESTART_NAME);
        std::cout << "> Initialization.name = " << restart_filename << std::endl;
    }
    p_phi_equil = config->get_qualified_as<double>("Initialization.p_phi_equil").value_or(DEFAULT_P_STOICH);
    std::cout << "> Initialization.p_phi_equil = " << p_phi_equil << std::endl;

    // Numerics
    n_particles = config->get_qualified_as<unsigned int>("Numerics.n_particles").value_or(DEFAULT_N_PARTICLES);
    std::cout << "> Numerics.n_particles = " << n_particles << std::endl;
    n_steps = config->get_qualified_as<unsigned int>("Numerics.n_steps").value_or(DEFAULT_N_STEPS);
    if (n_steps > 0)
        std::cout << "> Numerics.n_steps = " << n_steps << std::endl;
    t_stop = config->get_qualified_as<double>("Numerics.t_stop").value_or(DEFAULT_T_STOP);
    if (t_stop > 0.0)
        std::cout << "> Numerics.t_stop = " << t_stop << std::endl;
    dt_step = config->get_qualified_as<double>("Numerics.dt").value_or(DEFAULT_DT);
    if (dt_step > 0.0)
        std::cout << "> Numerics.dt = " << dt_step << std::endl;
    dt_sub_target = config->get_qualified_as<double>("Numerics.dt_sub").value_or(DEFAULT_DT_SUB);
    if (dt_sub_target > 0.0)
        std::cout << "> Numerics.dt_sub = " << dt_sub_target << std::endl;
    std::string convergence_metric_str =
        config->get_qualified_as<std::string>("Numerics.convergence_metric").value_or(DEFAULT_CONVERGENCE_METRIC);
    if (convergence_metric_str == "MEAN") {
        convergence_metric = MEAN;
    } else if (convergence_metric_str == "MEAN_VAR") {
        convergence_metric = MEAN_VAR;
    } else if (convergence_metric_str == "HIST") {
        convergence_metric = HIST;
    } else {
        std::cerr <<
          "Invalid convergence metric: " +
          convergence_metric_str +
          ". Must be MEAN, MEAN_VAR, or HIST." << std::endl;
        throw(0);
    }
    std::cout << "> Numerics.convergence_metric = " << convergenceMetricString(convergence_metric) << std::endl;

    n_stat = config->get_qualified_as<unsigned int>("Numerics.stats_window").value_or(DEFAULT_STATS_WINDOW);
    std::cout << "> Numerics.stats_window = " << n_stat << std::endl;
    rtol = config->get_qualified_as<double>("Numerics.rtol").value_or(DEFAULT_RTOL);
    if (rtol > 0.0)
        std::cout << "> Numerics.rtol = " << rtol << std::endl;
    min_steps_converge = config->get_qualified_as<unsigned int>("Numerics.min_steps_converge").value_or(DEFAULT_MIN_STEPS_CONVERGE);
    if (min_steps_converge > 0)
        std::cout << "> Numerics.min_steps_converge = " << min_steps_converge << std::endl;

    // Models
    std::string mixing_model_str = config->get_qualified_as<std::string>("Models.mixing_model").value_or(DEFAULT_MIX_MODEL);
    if (mixing_model_str == "NO_MIX") {
        mixing_model = NO_MIX;
    } else if (mixing_model_str == "FULL_MIX") {
        mixing_model = FULL_MIX;
    } else if (mixing_model_str == "CURL") {
        mixing_model = CURL;
    } else if (mixing_model_str == "MOD_CURL") {
        mixing_model = MOD_CURL;
    } else if (mixing_model_str == "IEM") {
        mixing_model = IEM;
    } else if (mixing_model_str == "EMST_1D") {
        mixing_model = EMST_1D;
    } else if (mixing_model_str == "EMST") {
        mixing_model = EMST;
    } else if (mixing_model_str == "KER_M") {
        mixing_model = KER_M;
    } else {
        std::cerr <<
          "Invalid mixing model: " +
          mixing_model_str +
          ". Must be NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, EMST_1D, EMST, or KER_M." << std::endl;
        throw(0);
    }
    std::cout << "> Models.mixing_model = " << mixingModelString(mixing_model) << std::endl;

    // Conditions
    std::string injection_mode_str =
        config->get_qualified_as<std::string>("Conditions.injection_mode").value_or(DEFAULT_INJECTION_MODE);
    if (injection_mode_str == "PREMIXED") {
        injection_mode = PREMIXED;
    } else if (injection_mode_str == "NONPREMIXED") {
        injection_mode = NONPREMIXED;
    } else {
        std::cerr <<
          "Invalid injection mode: " +
          injection_mode_str +
          ". Must be PREMIXED or NONPREMIXED." << std::endl;
        throw(0);
    }
    std::cout << "> Conditions.injection_mode = " << injectionModeString(injection_mode) << std::endl;
    pilot_flow = config->get_qualified_as<double>("Conditions.pilot_flow").value_or(DEFAULT_PILOT_FLOW);
    std::cout << "> Conditions.pilot_flow = " << pilot_flow << std::endl;
    P = config->get_qualified_as<double>("Conditions.pressure").value_or(DEFAULT_PRESSURE);
    std::cout << "> Conditions.pressure = " << P << std::endl;
    comp_fuel = config->get_qualified_as<std::string>("Conditions.comp_fuel").value_or(DEFAULT_COMP_FUEL);
    std::cout << "> Conditions.comp_fuel = " << comp_fuel << std::endl;
    comp_ox = config->get_qualified_as<std::string>("Conditions.comp_ox").value_or(DEFAULT_COMP_OX);
    std::cout << "> Conditions.comp_ox = " << comp_ox << std::endl;
    T_fuel = config->get_qualified_as<double>("Conditions.T_fuel").value_or(DEFAULT_T_FUEL);
    std::cout << "> Conditions.T_fuel = " << T_fuel << std::endl;
    T_ox = config->get_qualified_as<double>("Conditions.T_ox").value_or(DEFAULT_T_OX);
    std::cout << "> Conditions.T_ox = " << T_ox << std::endl;
    phi_global = config->get_qualified_as<double>("Conditions.phi_global").value_or(DEFAULT_PHI_GLOBAL);
    std::cout << "> Conditions.phi_global = " << phi_global << std::endl;
    std::string tau_res_mode_str = config->get_qualified_as<std::string>("Conditions.tau_res_mode").value_or(DEFAULT_TAU_RES_MODE);
    if (tau_res_mode_str == "EXP_MEAN") {
        tau_res_mode = EXP_MEAN;
    } else if (tau_res_mode_str == "DISTRIBUTION") {
        tau_res_mode = DISTRIBUTION;
    } else {
        std::cerr <<
          "Invalid tau_res mode: " +
          tau_res_mode_str +
          ". Must be EXP_MEAN or DISTRIBUTION." << std::endl;
        throw(0);
    }
    std::cout << "> Conditions.tau_res_Mode = " << tauResModeString(tau_res_mode) << std::endl;
    switch(tau_res_mode) {
        case EXP_MEAN: {
            tau_res_value = config->get_qualified_as<double>("Conditions.tau_res").value_or(DEFAULT_TAU_RES_VALUE);
            std::cout << "> Conditions.tau_res = " << tau_res_value << std::endl;
            break;
        }
        case DISTRIBUTION: {
            tau_res_hist_name = config->get_qualified_as<std::string>("Conditions.tau_res").value_or(DEFAULT_TAU_RES_HIST_NAME);
            std::cout << "> Conditions.tau_res = " << tau_res_hist_name << std::endl;
            break;
        }
        default: {
            break;
        }
    }
    tau_mix = config->get_qualified_as<double>("Conditions.tau_mix").value_or(DEFAULT_TAU_MIX);
    std::cout << "> Conditions.tau_mix = " << tau_mix << std::endl;

    // Output
    check_interval = config->get_qualified_as<unsigned int>("Output.check_interval").value_or(DEFAULT_CHECK_INTERVAL);
    std::cout << "> Output.check_interval = " << check_interval << std::endl;
    std::vector<std::string> check_variable_names_ = config->get_qualified_array_of<std::string>(
        "Output.check_variable_names").value_or(DEFAULT_CHECK_VARIABLE_NAMES);
    if (check_variable_names_.size() > 0) {
        check_variable_names.assign(check_variable_names_.begin(), check_variable_names_.end());
        std::cout << "> Output.check_variable_names: ";
        std::for_each(check_variable_names.begin(), check_variable_names.end(), [](std::string n){ std::cout << n << ", "; });
        std::cout << std::endl;
    }
    check_verbose = config->get_qualified_as<bool>("Output.check_verbose").value_or(DEFAULT_CHECK_VERBOSE);
    std::cout << "> Output.check_verbose = " << check_verbose << std::endl;
    write_raw_interval = config->get_qualified_as<unsigned int>("Output.write_raw_interval").value_or(DEFAULT_WRITE_INTERVAL);
    if (write_raw_interval < 0) {
        write_raw_interval = n_stat;
    }
    std::cout << "> Output.write_raw_interval = " << write_raw_interval << std::endl;
    if (write_raw_interval < n_stat) {
        std::cout << "WARNING: write_raw_interval < stats_window will result in duplicate data. " <<
            "Recommend write_raw_interval = stats_window" << std::endl;
    }
    write_stats_interval = config->get_qualified_as<unsigned int>("Output.write_stats_interval").value_or(DEFAULT_WRITE_INTERVAL);
    if (write_stats_interval < 0) write_stats_interval = check_interval;
    std::cout << "> Output.write_stats_interval = " << write_stats_interval << std::endl;

    std::cout << "--------------------------------------------------" << std::endl;
}

void PartiallyStirredReactor::initialize() {
    step = 0;
    t = 0.0;
    p_out = 0.0;
    p_mix = 0.0;
    rerror = std::numeric_limits<double>::infinity();

    // Create residence time histogram, if necessary
    if (tau_res_mode == DISTRIBUTION) {
        tau_res_hist.readHist(tau_res_hist_name);
        tau_res_hist.generatePDF();
        tau_res_hist.generateCDF();
    }

    // Set step sizes
    if (dt_step <= 0.0) {
        switch (tau_res_mode) {
            case EXP_MEAN: {
                dt_step = 0.1 * std::min(tau_res_value, tau_mix);
                break;
            }
            case DISTRIBUTION: {
                double tau_res_10 = tau_res_hist.percentileToValue(0.1);
                std::cout << "10th percentile residence time: " << tau_res_10 << std::endl;
                dt_step = 0.1 * std::min(tau_res_10, tau_mix);
                break;
            }
            // ^^ Capture the 10th percentile particle residence time
        }

        // EMST mixing models require a smaller time step
        if ((mixing_model == EMST_1D) || (mixing_model == EMST)) {
            dt_step /= 5;
        }

        std::cout << "Setting dt automatically: " << dt_step << std::endl;
    }
    if (dt_sub_target < 0.0) {
        dt_sub_target = 0.04 * tau_mix;
        std::cout << "Setting dt_sub automatically: " << dt_sub_target << std::endl;
    }

    // Create distributions and random engines for each thread
    for (int it = 0; it < omp_get_max_threads(); it++) {
        std::uniform_int_distribution<unsigned int> uni_int(0, n_particles-1);
        std::uniform_real_distribution<double> uni_real(0, 1);
        std::mt19937 eng((it + 1) * static_cast<uint64_t>(
            std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())));
        dists_uni_int.push_back(uni_int);
        dists_uni_real.push_back(uni_real);
        rand_engines.push_back(eng);
    }

    // Initialize ReactorNet for each thread
    std::cout << "Initializing ReactorNet for each thread..." << std::endl;
    solvec.resize(omp_get_max_threads());
    gasvec.resize(omp_get_max_threads());
    reactorvec.resize(omp_get_max_threads());
    rnetvec.resize(omp_get_max_threads());
#pragma omp parallel for
    for (int it = 0; it < omp_get_max_threads(); it++) {
        solvec[it] = Cantera::newSolution(mech_filename);
        gasvec[it] = solvec[it]->thermo();
        reactorvec[it] = new Cantera::IdealGasConstPressureReactor();
        reactorvec[it]->insert(solvec[it]);
        rnetvec[it] = new Cantera::ReactorNet();
        rnetvec[it]->addReactor(*reactorvec[it]);
    }
    std::cout << "Done initializing ReactorNets." << std::endl;
    
    n_species = gasvec[0]->nSpecies();
    n_state_variables = n_species + 1;
    pvec.resize(n_particles * n_stat);
    pvec_temp1.resize(n_particles);
    pvec_temp2.resize(n_particles);
    pvec_partemp.resize(omp_get_max_threads());
    xtemp1.resize(n_state_variables);
    xtemp2.resize(n_state_variables);
    xmean_old.resize(n_state_variables);
    xvar_old.resize(n_state_variables);
    Y_fuel.resize(n_species);
    Y_ox.resize(n_species);
    Y_mix.resize(n_species);

    // Get compositions
    gasvec[0]->setState_TPX(T_fuel, P, comp_fuel);
    gasvec[0]->getMassFractions(Y_fuel.data());
    h_fuel = gasvec[0]->enthalpy_mass();
    for (int k = 0; k < n_species; k++) {
        if (Y_fuel[k] > 1.0e-14) {
            i_fuel.push_back(k);
        }
    }
    gasvec[0]->setState_TPX(T_ox, P, comp_ox);
    gasvec[0]->getMassFractions(Y_ox.data());
    h_ox = gasvec[0]->enthalpy_mass();
    gasvec[0]->setEquivalenceRatio(phi_global, comp_fuel, comp_ox);
    gasvec[0]->getMassFractions(Y_mix.data());
    double Zeq = gasvec[0]->mixtureFraction(comp_fuel, comp_ox);
    h_mix = h_fuel * Zeq + h_ox * (1 - Zeq);

    // Compute equilibrium state (for initialization)
    gasvec[0]->setState_PY(P, Y_mix.data());
    gasvec[0]->setState_HP(h_mix, P);
    gasvec[0]->equilibrate("HP");
    double T_equil = gasvec[0]->temperature();
    double h_equil = gasvec[0]->enthalpy_mass();
    std::vector<double> Y_equil(n_species);
    gasvec[0]->getMassFractions(Y_equil.data());

    // Initialize particles
    std::cout << "Initializing particles..." << std::endl;
    int np_phi_equil = std::round(p_phi_equil * n_particles);

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for
        for (int ip = 0; ip < n_particles; ip++) {
            pvec[ip].setP(&P);
            pvec[ip].setID(id_iterator++);
            pvec[ip].setnSpecies(n_species);
            pvec[ip].setMass(1.0); // TODO: check this (fuel particles lighter?)

            if (ip < np_phi_equil) {
                pvec[ip].seth(h_equil);
                pvec[ip].setY(Y_equil.data());
            } else {
                double Z = (double)(ip - np_phi_equil) / (double)(n_particles - np_phi_equil);
                gasvec[tid]->setMixtureFraction(Z, comp_fuel, comp_ox);
                gasvec[tid]->setState_HP(h_mix, P);
                gasvec[tid]->equilibrate("HP");
                pvec[ip].seth(gasvec[tid]->enthalpy_mass());
                pvec[ip].setY(gasvec[tid]->massFractions());
            }

            switch (tau_res_mode) {
                case EXP_MEAN: {
                    pvec[ip].setAge(0.0);
                    pvec[ip].setTauRes(0.0); // TODO - figure out how to compute this in this mode
                    break;
                }
                case DISTRIBUTION: {
                    pvec[ip].setAge(0.0);
                    pvec[ip].setTauRes(tau_res_hist.rand(dists_uni_real[tid], rand_engines[tid]));
                    break;
                }
            }
        }
    }

    // Read restart
    if (restart) {
        std::cout << "Reading restart..." << std::endl;
        readRestart();
        
        std::cout << "Done reading restart." << std::endl;
    }

    // Copy state to all old particles to initialize stats properly
#pragma omp parallel for
    for (int ip = 0; ip < n_particles; ip++) {
        for (int is = 1; is < n_stat; is++) {
            pvec[ip + is*n_particles] = pvec[ip];
        }
    }
    std::cout << "Done initializing particles." << std::endl;

    // Initialize injectors
    injvec.push_back(Injector(0, n_species));
    injvec[0].seth(h_equil);
    injvec[0].setY(Y_equil.data());
    injvec[0].setFlow(pilot_flow);

    switch (injection_mode) {
        case PREMIXED: {
            injvec.push_back(Injector(1, n_species));

            // Premix injector
            injvec[1].seth(h_mix);
            injvec[1].setY(Y_mix.data());
            injvec[1].setFlow(1.0 - pilot_flow);
            break;
        }
        case NONPREMIXED: {
            for (int iinj = 1; iinj < 3; iinj++) {
                injvec.push_back(Injector(iinj, n_species));
            }

            // Fuel injector
            injvec[1].seth(h_fuel);
            injvec[1].setY(Y_fuel.data());
            injvec[1].setFlow(Zeq * (1.0-pilot_flow));

            // Oxidizer injector
            injvec[2].seth(h_ox);
            injvec[2].setY(Y_ox.data());
            injvec[2].setFlow((1-Zeq) * (1.0-pilot_flow));
            break;
        }
        default: {
            throw Cantera::CanteraError("PartiallyStirredReactor::subStepMix",
                                        "Invalid mixing model.");
        }
    }

    for (Injector& inj : injvec) {
        inj.print();
    }

    // Initialize variable functions
    for (int iv = 0; iv < n_state_variables; iv++) {
        variable_functions.push_back(
            [this, iv](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
                return pvec[ip].state(iv); });
    }

    // Initialize auxiliary variables
    aux_variable_names.push_back("id");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getID(); });
    n_aux_variables++;

    aux_variable_names.push_back("age");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getAge(); });
    n_aux_variables++;

    aux_variable_names.push_back("tau_res");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getTauRes(); });
    n_aux_variables++;

    aux_variable_names.push_back("mass");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getMass(); });
    n_aux_variables++;

    // Initialize derived variables
    derived_variable_names.push_back("T");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].T(gas); });
    n_derived_variables++;

    derived_variable_names.push_back("Z");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].Z(gas, comp_fuel, comp_ox); });
    n_derived_variables++;

    // Set check variables
    if (check_variable_names.size() > 0) {
        for (auto& name : check_variable_names) {
            iv_check.push_back(variableIndex(name));
        }
    } else {
        for (int iv = 0; iv < nVariables(); iv++) {
            iv_check.push_back(iv);
        }
    }

    // Initialize statistics
    meanState(&xmean_old, true, true);
    varianceState(&xvar_old, true, true);

    // Initialize output files
    writeRawHeaders();
    writeStatsHeaders();
}

void PartiallyStirredReactor::readRestart() {
    std::string line, word;
    std::ifstream file(restart_filename);
    std::vector<int> iv_read;
    if (file.is_open()) {

        std::stringstream str(line);

        // Read header
        std::getline(file, line);
        while (std::getline(str, word, ',')) {
            iv_read.push_back(variableIndex(word));
        }

        // Read content
        for (int ip = 0; ip < n_particles; ip++) {
            if (std::getline(file, line)) {
                for (int iiv = 0; iiv < iv_read.size(); iiv++) {
                    if (iv_read[iiv] >= n_state_variables + n_aux_variables) continue;
                    std::getline(str, word, ',');
                    pvec[ip].stateAux(iv_read[iiv]) = std::stod(word);
                }
            } else {
                std::cout << "WARNING: Only " << ip+1 << " particle(s) specified in restart file. Repeating read." << std::endl;
                file.clear(); // Reset eof and fail flags
                file.seekg(0); // Go to beginning of file
                std::getline(file, line); // Advance by 1 line
                ip--; // Decrement ip so this particle is not skipp
            }
        }

        // Warn about extra particles
        if (std::getline(file, line)) {
            std::cout << "WARNING: More particles than necessary specified in restart file. " << 
                "Only reading first " << n_particles << "particle(s)." << std::endl;
        }
    } else {
        throw std::runtime_error("PartiallyStirredReactor::readRestart - could not open " + restart_filename + ".");
    }
}

void PartiallyStirredReactor::run() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Begin time stepping..." << std::endl;
    run_done = false;
    particles_injected = false;
    while (!run_done) {
        takeStep();
    }
    std::cout << "Done." << std::endl;
}

void PartiallyStirredReactor::print() {
    
}

void PartiallyStirredReactor::checkVariable(int iv) {
    double meanval = mean(iv, true, true);
    std::cout <<
        std::left << std::setw(COL_WIDTH) << "> " + variableName(iv) <<
        std::left << std::setw(COL_WIDTH) << min(iv, true) <<
        std::left << std::setw(COL_WIDTH) << meanval <<
        std::left << std::setw(COL_WIDTH) << variance(iv, meanval, true) <<
        std::left << std::setw(COL_WIDTH) << max(iv, true) << std::endl;
}

void PartiallyStirredReactor::check(bool force) {
    if (!force) {
        if (check_interval == 0) return;
        else if ((step % check_interval) != 0) return;
    }
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Starting step: " << step << ", t = " << t << std::endl;
    std::cout << std::endl;
    std::cout << "--- Stats for last " << n_stat <<
        " steps as of beginning of step " << step << " ---" << std::endl;
    std::cout <<
        std::left << std::setw(COL_WIDTH) << "> Name" <<
        std::left << std::setw(COL_WIDTH) << "Min" <<
        std::left << std::setw(COL_WIDTH) << "Favre Mean" <<
        std::left << std::setw(COL_WIDTH) << "Variance" <<
        std::left << std::setw(COL_WIDTH) << "Max" << std::endl;
    if (check_verbose) {
        for (int iv = 0; iv < nVariables(); iv++) {
            checkVariable(iv);
        }
    } else {
        for (auto& iv : iv_check) {
            checkVariable(iv);
        }
    }

    std::cout << std::endl;
    std::cout << "> n_recycled = " << n_recycled << " last step, " <<
        n_recycled_check << " since last check" << std::endl;
    std::cout << "> rerror = " << rerror << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    n_recycled_check = 0;
}

std::string PartiallyStirredReactor::variableName(int iv) {
    if (iv == c_offset_h) {
        return "h";
    } else if (iv < n_state_variables) {
        std::string prefix = "Y_";
        return prefix + gasvec[0]->speciesName(iv-c_offset_Y);
    } else if (iv < n_state_variables + n_aux_variables) {
        return aux_variable_names[iv - n_state_variables];
    } else {
        return derived_variable_names[iv - n_state_variables - n_aux_variables];
    }
}

int PartiallyStirredReactor::variableIndex(std::string name) {
    for (int iv = 0; iv < nVariables(); iv++) {
        if (name == variableName(iv)) {
            return iv;
        }
    }
    throw std::runtime_error("PartiallyStirredReactor::variableIndex - variable " + name + " not found.");
}

void PartiallyStirredReactor::takeStep() {
    // Print status
    check();

    // Copy state to other stats region
    copyState();

    // Adjust time step
    calcDt();

    // Take substeps
    // for (Particle& p : pvec) p.print(1.0e-14, gasvec[0]);
    subStepInflow(dt_step);
    // for (Particle& p : pvec) p.print(1.0e-14, gasvec[0]);
    for (int isub = 0; isub < n_sub; isub++) {
        subStepMix(dt_sub);
        // for (Particle& p : pvec) p.print(1.0e-14, gasvec[0]);
        subStepReact(dt_sub);
        // for (Particle& p : pvec) p.print(1.0e-14, gasvec[0]);
    }
    // for (Particle& p : pvec) p.print(1.0e-14, gasvec[0]);

    // Increment counters
    incrementAge();
    step++;
    t += dt_step;

    // Calculate convergence metric
    calcConvergence();

    // Determine if run is complete
    run_done = runDone();

    // Write data (force if run is complete)
    writeRaw(run_done);
    writeStats(run_done);
}

void PartiallyStirredReactor::calcDt() {
    if ((t_stop > 0) && ((t + dt_step) > t_stop)) {
        dt_step = t_stop - t;
    }
    n_sub = 1 + std::round(dt_step / dt_sub_target);
    // n_sub = 1; // DEBUG
    dt_sub = dt_step / n_sub;
}

void PartiallyStirredReactor::subStepInflow(double dt) {
    int n_recycled_ = 0;
    int n_recycled_check_ = n_recycled_check;
    switch (tau_res_mode) {
        case EXP_MEAN: {
            // Random choice formulation
            p_out += n_particles * dt / tau_res_value; // Fractional particle count to recycle
            int np_out = std::round(p_out); // Round to integer
            p_out -= np_out; // Hold on to remainder for next step

            std::vector<int> ip_vec(np_out, -1);
            bool repeat;
            for (int i = 0; i < np_out; i++) {

                // Generate random particle index to recycle
                unsigned int ip = dists_uni_int[0](rand_engines[0]);

                // Check if repeated
                repeat = false;
                for (int j = 0; j < i; j++) {
                    if (ip_vec[j] == ip) {
                        repeat = true;
                        break;
                    }
                }

                if (repeat) {
                    i--; // If repeated, try again
                } else {
                    ip_vec[i] = ip; // Otherwise, set
                }
            }

#pragma omp parallel reduction(+:n_recycled_,n_recycled_check_)
            {
                int tid = omp_get_thread_num();
#pragma omp for
                // Iterate over particles and recycle
                for (int ip_out = 0; ip_out < np_out; ip_out++) {
                    double p_inj = dists_uni_real[tid](rand_engines[tid]); // Probability: which injector to use for new particle
                    recycleParticle(ip_vec[ip_out], p_inj, tid);
                    n_recycled_++;
                    n_recycled_check_++;
                }
            }
            break;
        }
        case DISTRIBUTION: {
            // Life expectancy formulation
#pragma omp parallel reduction(+:n_recycled_,n_recycled_check_)
            {
                int tid = omp_get_thread_num();
#pragma omp for
                // Iterate over particles and recycle
                for (int ip = 0; ip < n_particles; ip++) {
                    if (pvec[ip].tooOld()) {
                        double p_inj = dists_uni_real[tid](rand_engines[tid]); // Probability: which injector to use for new particle
                        recycleParticle(ip, p_inj, tid);
                        n_recycled_++;
                        n_recycled_check_++;
                    }
                }
            }
            break;
        }
    }
    n_recycled = n_recycled_;
    n_recycled_check = n_recycled_check_;
    particles_injected = particles_injected || (n_recycled > 0);
}

void PartiallyStirredReactor::subStepMix(double dt) {
    switch(mixing_model) {
        case NO_MIX: {
            // Do nothing
            break;
        }
        case FULL_MIX: {
            meanState(&xtemp1, false, true);
            // Set all particles to Favre mean state
#pragma omp parallel for
            for (int ip = 0; ip < n_particles; ip++) {
                pvec[ip].setState(xtemp1.data());
            }
            break;
        }
        case CURL: {
            // Compute how many pairs to mix
            p_mix += n_particles * dt / tau_mix;
            int np_mix = std::round(p_mix);
            p_mix -= np_mix;
#pragma omp parallel
            {
                int tid = omp_get_thread_num();
#pragma omp for
                // Iterate over pairs and mix
                for (int ipair = 0; ipair < np_mix; ipair++) {
                    // TODO: particle weights are not used
                    unsigned int ip1 = dists_uni_int[tid](rand_engines[tid]);
                    unsigned int ip2 = dists_uni_int[tid](rand_engines[tid]);
                    double a = dists_uni_real[tid](rand_engines[tid]);
                    pvec_partemp[tid] = (pvec[ip1].getMass() * pvec[ip1] + pvec[ip2].getMass() * pvec[ip1]) /
                        (pvec[ip1].getMass() + pvec[ip2].getMass());
                    pvec[ip1] = pvec_partemp[tid];
                    pvec[ip2] = pvec_partemp[tid];
                }
            }
            break;
        }
        case MOD_CURL: {
            // Compute how many pairs to mix
            p_mix += n_particles * dt / tau_mix;
            int np_mix = std::round(p_mix);
            p_mix -= np_mix;
#pragma omp parallel
            {
                int tid = omp_get_thread_num();
#pragma omp for
                // Iterate over pairs and mix
                for (int ipair = 0; ipair < np_mix; ipair++) {
                    // TODO: Particle weights are not used
                    unsigned int ip1 = dists_uni_int[tid](rand_engines[tid]);
                    unsigned int ip2 = dists_uni_int[tid](rand_engines[tid]);
                    double a = dists_uni_real[tid](rand_engines[tid]);
                    pvec_partemp[tid] = (pvec[ip1].getMass() * pvec[ip1] + pvec[ip2].getMass() * pvec[ip1]) /
                        (pvec[ip1].getMass() + pvec[ip2].getMass());
                    pvec[ip1] += a * (pvec_partemp[tid] - pvec[ip1]);
                    pvec[ip2] += a * (pvec_partemp[tid] - pvec[ip2]);
                }
            }
            break;
        }
        case IEM: {
            // meanState(&xtemp1, false, true);
            meanState(&xtemp1, false, false); // DEBUG
            Particle pmean;
            pmean.setnSpecies(n_species);
            pmean.setState(xtemp1.data());
#pragma omp parallel for
            for (int ip = 0; ip < n_particles; ip++) {
                pvec[ip] += -(1.0/2.0) * (dt/tau_mix) * (pvec[ip] - pmean);
            }
            break;
        }
        case EMST_1D: {
            int ip_m, ip_n;
            double Wi, Wv;
            double dt_save = dt;
            double dt_temp = dt;
            double alpha;
            std::vector<double> w(n_particles, 0.0);
            std::vector<double> Bv(n_particles, 0.0);
            std::vector<double> A(n_state_variables, 0.0);
            std::vector<double> B(n_state_variables, 0.0);
            std::vector<double> C(n_state_variables, 0.0);

            // Get particle weights
            double M_tot = 0.0;
#pragma omp parallel for reduction(+:M_tot)
            for (int ip = 0; ip < n_particles; ip++) {
                w[ip] = pvec[ip].getMass();
                M_tot += w[ip];
            }

#pragma omp parallel for
            for (int ip = 0; ip < n_particles; ip++) {
                w[ip] /= M_tot;
            }

            while (dt_save > 0) {

                // Get particle mixture fractions and masses
                std::vector<double> Z(n_particles, 0.0);
#pragma omp parallel for
                for (int ip = 0; ip < n_particles; ip++) {
                    Z[ip] = pvec[ip].Z(gasvec[omp_get_thread_num()], comp_fuel, comp_ox);
                }

                // Get mixture fraction sorting index
                std::vector<std::size_t> iZ_sorted = sort_indices(Z);

                // Compute variance of current state
                varianceState(&xtemp1, false, false);

                // Compute weights
                Wi = 0;
                for (int ip = 0; ip < n_particles; ip++) {
                    Wi += w[iZ_sorted[ip]];
                    Wv = std::min(Wi, 1 - Wi);
                    Bv[ip] = 2 * Wv;
                }

                // Initialize temporary particle vector
#pragma omp parallel for
                for (int ip = 0; ip < n_particles; ip++) {
                    pvec_temp1[ip] = pvec[ip];
                    for (int iv = 0; iv < n_state_variables; iv++) {
                        pvec_temp1[ip].state(iv) = 0.0;
                    }
                }

                // Compute mixing dpvec
// #pragma omp parallel for reduction(vec_particle_plus:pvec_temp1)
                for (int ip = 0; ip < n_particles - 1; ip++) {
                    ip_m = iZ_sorted[ip];
                    ip_n = iZ_sorted[ip+1];
                    pvec_temp1[ip_m] += - Bv[ip] * (pvec[ip_m] - pvec[ip_n]) / w[ip_m];
                    pvec_temp1[ip_n] += - Bv[ip] * (pvec[ip_n] - pvec[ip_m]) / w[ip_n];
                }

                // Determine maximum allowable time step
                for (int iv = 0; iv < n_state_variables; iv++) {
                    double A_ = 0.0;
                    double B_ = 0.0;
                    double C_ = 0.0;
#pragma omp parallel for reduction(+:A_,B_,C_)
                    for (int ip = 0; ip < n_particles; ip++) {
                        A_ += pvec_temp1[ip].state(iv) * pvec_temp1[ip].state(iv);
                        B_ += pvec_temp1[ip].state(iv) * pvec[ip].state(iv);
                        C_ += (1.0 / tau_mix) * xtemp1[iv];
                    }
                    A[iv] = A_ / n_particles;
                    B[iv] = 2.0 * B_ / n_particles;
                    C[iv] = C_;
                    dt_temp = std::min(dt_temp, (B[iv] * B[iv]) / (4.0 * A[iv] * C[iv]));
                }

                alpha = std::numeric_limits<double>::infinity();
#pragma omp parallel for reduction(min:alpha)
                for (int iv = 0; iv < n_state_variables; iv++) {
                    alpha = std::min(alpha, -B[iv] / (2 * A[iv] * dt_temp));
                }

                dt_temp = std::min(dt_temp, dt_save);

                // Root finding - search for optimal mixing ratio
                for (int ii = 0; ii < 4; ii++) {
#pragma omp parallel for
                    for (int ip = 0; ip < n_particles; ip++) {
                        pvec_temp2[ip] = pvec[ip] + pvec_temp1[ip] * alpha * dt_temp;
                    }
                    varianceState(&pvec_temp2, &xtemp2, false);
                    double var_decay = 0.0;
                    for (int iv = 0; iv < n_state_variables; iv++) {
                        double var_decay_iv = xtemp2[iv] / xtemp1[iv];
                        // Check NaN, arises when variance is 0
                        if (var_decay_iv == var_decay_iv) {
                            var_decay += var_decay_iv;
                        }

                        // DEBUG
                        // std::cout << "iv: " << iv <<
                        //     " xtemp1[iv]: " << xtemp1[iv] <<
                        //     " xtemp2[iv]: " << xtemp2[iv] <<
                        //     " var_decay: " << var_decay << std::endl;
                        // DEBUG
                    }
                    var_decay = 1.0 - (var_decay / n_state_variables);
                    double var_ratio = var_decay / (1.0 - std::exp(-dt_temp / tau_mix));
                    alpha /= var_ratio;
                }

                // Update state
#pragma omp parallel for
                for (int ip = 0; ip < n_particles; ip++) {
                    pvec[ip] += pvec_temp1[ip] * alpha * dt_temp;
                }

                dt_save -= dt_temp;
            }
            break;
        }
        case EMST: {
            throw Cantera::NotImplementedError("PartiallyStirredReactor::subStepMix",
                                               "EMST not implemented.");
            break;
        }
        case KER_M: {
            throw Cantera::NotImplementedError("PartiallyStirredReactor::subStepMix",
                                               "KER_M not implemented.");
            break;
        }
        default: {
            throw Cantera::CanteraError("PartiallyStirredReactor::subStepMix",
                                        "Invalid mixing model.");
        }
    }
}

void PartiallyStirredReactor::subStepReact(double dt) {
#pragma omp parallel for
    for (int ip = 0; ip < n_particles; ip++) {
        // Use the current thread's reactor for calculation
        pvec[ip].react(rnetvec[omp_get_thread_num()], dt);
    }
}

void PartiallyStirredReactor::incrementAge() {
#pragma omp parallel for
    for (int ip = 0; ip < n_particles; ip++) {
        pvec[ip].getAge() += dt_step;
    }
}

void PartiallyStirredReactor::recycleParticle(unsigned int ip, double p_inj, int tid) {
    int iinj;
    double flow_sum = 0.0;
    for (int i = 0; i < injvec.size(); i++) {
        flow_sum += injvec[i].getFlow();
        if (p_inj <= flow_sum) {
            iinj = i;
            break;
        }
    }
    pvec[ip].setID(id_iterator++);
    pvec[ip].setAge(0.0);
    switch (tau_res_mode) {
        case EXP_MEAN: {
            pvec[ip].setTauRes(0.0);
            break;
        }
        case DISTRIBUTION: {
            pvec[ip].setTauRes(tau_res_hist.rand(dists_uni_real[tid], rand_engines[tid]));
            break;
        }
    }
    pvec[ip].seth(injvec[iinj].h());
    pvec[ip].setY(injvec[iinj].Y().data());
    pvec[ip].getnRecycles()++;
}

void PartiallyStirredReactor::calcConvergence() {
    switch (convergence_metric) {
        case MEAN: {
            meanState(&xtemp1, true, true);
            rerror = 0.0;
            for (int iv = 0; iv < n_state_variables; iv++) {
                rerror += std::pow((xtemp1[iv] - xmean_old[iv]), 2.0);
            }
            rerror = std::sqrt(rerror);
            break;
        }
        case MEAN_VAR: {
            meanState(&xtemp1, true, true);
            varianceState(&xtemp2, &xtemp1, true);
            rerror = 0.0;
            for (int iv = 0; iv < n_state_variables; iv++) {
                rerror += std::pow((xtemp1[iv] - xmean_old[iv]), 2.0);
                rerror += std::pow((xtemp2[iv] - xvar_old[iv]), 2.0);
            }
            rerror = std::sqrt(rerror);
            break;
        }
        case HIST: {
            throw Cantera::NotImplementedError("PartiallyStirredReactor::calcConvergence",
                                               "MEAN_VAR not implemented.");
            break;
        }
        default: {
            throw Cantera::CanteraError("PartiallyStirredReactor::calcConvergence",
                                        "Invalid mixing model.");
        }
    }
}

bool PartiallyStirredReactor::runDone() {
    if ((n_steps > 0) && (step >= n_steps)) {
        std::cout << "Reached termination condition: step >= n_steps at time t = " << t << std::endl;
        return true;
    }
    if ((t_stop > 0.0) && (t >= t_stop)) {
        std::cout << "Reached termination condition: t >= t_stop at step " << step << std::endl;
        return true;
    }
    if ((rerror <= rtol) && (step >= min_steps_converge) && (particles_injected)) {
        std::cout << "Reached termination condition: rerror (" <<
            rerror << ") <= rtol (" << rtol << ") at step " << step << std::endl;
        return true;
    }
    double fmeanT = mean(variableIndex("T"), true, true);
    if (fmeanT < T_EXTINCT) {
        std::cout << "Reached termination condition: fmean(T) (" <<
            fmeanT << ") <= T_EXTINCT (" << T_EXTINCT << ") at step " << step << std::endl;
        return true;
    }
    return false;
}

void PartiallyStirredReactor::copyState() {
    if (n_stat == 1) return;
#pragma omp parallel for
    for (int ip = 0; ip < n_particles; ip++) {
        pvec[ip + i_stat*n_particles] = pvec[ip];
    }
    if (i_stat >= n_stat-1) {
        i_stat = 1;
    } else {
        i_stat++;
    }
}

void PartiallyStirredReactor::writeRawHeaders() {
    std::ofstream file;
    file.open(RAW_NAME + RAW_EXT);
    for (int iv = 0; iv < nVariables(); iv++) {
        file << variableName(iv);
        if (iv < nVariables()-1) {
            file << ",";
        } else {
            file << std::endl;
        }
    }
    file.close();
}

void PartiallyStirredReactor::writeRaw(bool force) {
    if ((!force) && (step % write_raw_interval != 0)) return;

    std::ofstream file;
    file.open(RAW_NAME + RAW_EXT, std::ios_base::app);
    for (int ip = 0; ip < n_particles * n_stat; ip++) {
        for (int iv = 0; iv < nVariables(); iv++) {
            file << std::setprecision(WRITE_PRECISION) << variable_functions[iv](gasvec[0], ip);
            if (iv < nVariables()-1) {
                file << ",";
            } else {
                file << std::endl;
            }
        }
    }
    file.close();
}

std::string PartiallyStirredReactor::statsPath(std::string name) {
    return STATS_DIR + "/" + STATS_PREFIX + name + STATS_EXT;
}

void PartiallyStirredReactor::writeStatsHeaders() {
    fs::create_directories(STATS_DIR);

    std::ofstream file_min;
    std::ofstream file_max;
    std::ofstream file_fmean;
    std::ofstream file_variance;

    file_min.open(statsPath("min"));
    file_max.open(statsPath("max"));
    file_fmean.open(statsPath("fmean"));
    file_variance.open(statsPath("variance"));

    file_min << "step,time,";
    file_max << "step,time,";
    file_fmean << "step,time,";
    file_variance << "step,time,";

    for (int iv = 0; iv < nVariables(); iv++) {
        file_min << variableName(iv);
        file_max << variableName(iv);
        file_fmean << variableName(iv);
        file_variance << variableName(iv);
        if (iv < nVariables()-1) {
            file_min << ",";
            file_max << ",";
            file_fmean << ",";
            file_variance << ",";
        } else {
            file_min << std::endl;
            file_max << std::endl;
            file_fmean << std::endl;
            file_variance << std::endl;
        }
    }

    file_min.close();
    file_max.close();
    file_fmean.close();
    file_variance.close();
}

void PartiallyStirredReactor::writeStats(bool force) {
    if ((!force) && (step % write_stats_interval != 0)) return;

    std::ofstream file_min;
    std::ofstream file_max;
    std::ofstream file_fmean;
    std::ofstream file_variance;

    file_min.open(statsPath("min"), std::ios_base::app);
    file_max.open(statsPath("max"), std::ios_base::app);
    file_fmean.open(statsPath("fmean"), std::ios_base::app);
    file_variance.open(statsPath("variance"), std::ios_base::app);

    file_min << step << "," << t << ",";
    file_max << step << "," << t << ",";
    file_fmean << step << "," << t << ",";
    file_variance << step << "," << t << ",";

    for (int iv = 0; iv < nVariables(); iv++) {
        double meanval = mean(iv, true, true);
        file_min << std::setprecision(WRITE_PRECISION) << min(iv, true);
        file_max << std::setprecision(WRITE_PRECISION) << max(iv, true);
        file_fmean << std::setprecision(WRITE_PRECISION) << meanval;
        file_variance << std::setprecision(WRITE_PRECISION) << variance(iv, meanval, true);

        if (iv < nVariables()-1) {
            file_min << ",";
            file_max << ",";
            file_fmean << ",";
            file_variance << ",";
        } else {
            file_min << std::endl;
            file_max << std::endl;
            file_fmean << std::endl;
            file_variance << std::endl;
        }
    }

    file_min.close();
    file_max.close();
    file_fmean.close();
    file_variance.close();
}

double PartiallyStirredReactor::min(int iv, bool all) {
    return min(variable_functions[iv], all);
}

double PartiallyStirredReactor::max(int iv, bool all) {
    return max(variable_functions[iv], all);
}

double PartiallyStirredReactor::mean(int iv, bool all, bool favre) {
    return mean(variable_functions[iv], all, favre);
}

double PartiallyStirredReactor::variance(int iv, bool all, bool favre) {
    double meanval = mean(variable_functions[iv], all, favre);
    return variance(variable_functions[iv], meanval, all);
}

double PartiallyStirredReactor::variance(int iv, double meanval, bool all) {
    return variance(variable_functions[iv], meanval, all);
}

double PartiallyStirredReactor::min(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;
    double minval = std::numeric_limits<double>::infinity();
#pragma omp parallel for reduction(min:minval)
    for (int ip = 0; ip < ip_stop; ip++) {
        minval = std::min(minval, xfunc(gasvec[omp_get_thread_num()], ip));
    }
    return minval;
}

double PartiallyStirredReactor::max(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;
    double maxval = -std::numeric_limits<double>::infinity();
#pragma omp parallel for reduction(max:maxval)
    for (int ip = 0; ip < ip_stop; ip++) {
        maxval = std::max(maxval, xfunc(gasvec[omp_get_thread_num()], ip));
    }
    return maxval;
}

double PartiallyStirredReactor::mean(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all, bool favre) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;
    double rhosum = 0.0;
    double xsum = 0.0;
#pragma omp parallel for reduction(+:rhosum,xsum)
    for (int ip = 0; ip < ip_stop; ip++) {
        double rho = 0.0;
        if (favre) {
            rho = pvec[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        xsum += rho * xfunc(gasvec[omp_get_thread_num()], ip);
    }
    return xsum / rhosum;
}

double PartiallyStirredReactor::variance(std::function<double(std::shared_ptr<Cantera::ThermoPhase>, int)> xfunc, bool all, bool favre) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;
    double meanval = mean(xfunc, favre);
    double varsum = 0.0;
#pragma omp parallel for reduction(+:varsum)
    for (int ip = 0; ip < ip_stop; ip++) {
        varsum += std::pow((xfunc(gasvec[omp_get_thread_num()], ip) - meanval), 2.0);
    }
    return varsum / n_particles;
}

void PartiallyStirredReactor::minState(std::vector<double>* minvec, bool all) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;

    // Min across particles
    std::vector<double> minvec_temp(minvec->size(), std::numeric_limits<double>::infinity());
#pragma omp parallel for reduction(vec_double_min:minvec_temp)
    for (int ip = 0; ip < ip_stop; ip++) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            minvec_temp[iv] = std::min(minvec_temp[iv], pvec[ip].state(iv));
        }
    }

    // Write to minvec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*minvec)[iv] = minvec_temp[iv];
    }
}

void PartiallyStirredReactor::maxState(std::vector<double>* maxvec, bool all) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;

    // Max across particles
    std::vector<double> maxvec_temp(maxvec->size(), -std::numeric_limits<double>::infinity());
#pragma omp parallel for reduction(vec_double_max:maxvec_temp)
    for (int ip = 0; ip < ip_stop; ip++) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            maxvec_temp[iv] = std::max(maxvec_temp[iv], pvec[ip].state(iv));
        }
    }

    // Write to maxvec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*maxvec)[iv] = maxvec_temp[iv];
    }
}

void PartiallyStirredReactor::meanState(std::vector<double>* xmeanvec, bool all, bool favre) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;

    // Sum across particles
    std::vector<double> xsumvec(xmeanvec->size(), 0.0);
    double rhosum = 0.0;
#pragma omp parallel for reduction(+:rhosum) reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < ip_stop; ip++) {
        double rho;
        if (favre) {
            rho = pvec[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        for (int iv = 0; iv < n_state_variables; iv++) {
            xsumvec[iv] += rho * pvec[ip].state(iv);
        }
    }

    // Divide by particle count and write to mean vec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*xmeanvec)[iv] = xsumvec[iv] / rhosum;
    }
}

void PartiallyStirredReactor::meanState(std::vector<Particle>* pvec_, std::vector<double>* xmeanvec, bool favre) {
    // Sum across particles
    std::vector<double> xsumvec(xmeanvec->size(), 0.0);
    double rhosum = 0.0;
#pragma omp parallel for reduction(+:rhosum) reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < n_particles; ip++) {
        double rho;
        if (favre) {
            rho = (*pvec_)[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        for (int iv = 0; iv < n_state_variables; iv++) {
            xsumvec[iv] += rho * (*pvec_)[ip].state(iv);
        }
    }

    // Divide by particle count and write to mean vec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*xmeanvec)[iv] = xsumvec[iv] / rhosum;
    }
}

void PartiallyStirredReactor::varianceState(std::vector<double>* xvarvec, bool all, bool favre) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;

    // Sum across particles
    std::vector<double> xmeanvec(xvarvec->size(), 0.0);
    meanState(&xmeanvec, favre);
    std::vector<double> xsumvec(xvarvec->size(), 0.0);
#pragma omp parallel for reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < ip_stop; ip++) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xsumvec[iv] += std::pow((pvec[ip].state(iv) - xmeanvec[iv]), 2.0);
        }
    }

    // Divide by particle count and write to variance vec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*xvarvec)[iv] = xsumvec[iv] / n_particles;
    }
}

void PartiallyStirredReactor::varianceState(std::vector<Particle>* pvec_, std::vector<double>* xvarvec, bool favre) {

    // Sum across particles
    std::vector<double> xmeanvec(xvarvec->size(), 0.0);
    meanState(pvec_, &xmeanvec, favre);
    std::vector<double> xsumvec(xvarvec->size(), 0.0);
#pragma omp parallel for reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < n_particles; ip++) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xsumvec[iv] += std::pow(((*pvec_)[ip].state(iv) - xmeanvec[iv]), 2.0);
        }
    }

    // Divide by particle count and write to variance vec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*xvarvec)[iv] = xsumvec[iv] / n_particles;
    }
}

void PartiallyStirredReactor::varianceState(std::vector<double>* xvarvec, std::vector<double>* xmeanvec, bool all) {
    int ip_stop = (all) ? n_particles * n_stat : n_particles;

    // Sum across particles
    std::vector<double> xsumvec(xvarvec->size(), 0.0);
#pragma omp parallel for reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < ip_stop; ip++) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xsumvec[iv] += std::pow((pvec[ip].state(iv) - (*xmeanvec)[iv]), 2.0);
        }
    }

    // Divide by particle count and write to variance vec
    for (int iv = 0; iv < n_state_variables; iv++) {
        (*xvarvec)[iv] = xsumvec[iv] / n_particles;
    }
}

PartiallyStirredReactor::~PartiallyStirredReactor() {
}
