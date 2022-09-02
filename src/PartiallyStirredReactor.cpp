#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include "../cpptoml/include/cpptoml.h"

#include "common.h"
#include "PartiallyStirredReactor.h"

PartiallyStirredReactor::PartiallyStirredReactor(const std::string& input_filename_) :
    input_filename(input_filename_),
    step(0), t(0.0), p_out(0.0), i_stat(1), n_particles(0), n_state_variables(0), n_derived_variables(0), n_species(0)
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
    } else if (mixing_model_str == "EMST") {
        mixing_model = EMST;
    } else {
        std::cerr <<
          "Invalid mixing model: " +
          mixing_model_str +
          ". Must be NO_MIX, FULL_MIX, CURL, MOD_CURL, IEM, or EMST." << std::endl;
        throw(0);
    }
    std::cout << "> Models.mixing_model = " << mixingModelString(mixing_model) << std::endl;

    // Conditions
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
    T_init = config->get_qualified_as<double>("Conditions.T_init").value_or(DEFAULT_T_INIT);
    std::cout << "> Conditions.T_init = " << T_init << std::endl;
    phi_global = config->get_qualified_as<double>("Conditions.phi_global").value_or(DEFAULT_PHI_GLOBAL);
    std::cout << "> Conditions.phi_global = " << phi_global << std::endl;
    std::string tau_res_mode_str = config->get_qualified_as<std::string>("Conditions.tau_res_mode").value_or(DEFAULT_TAU_RES_MODE);
    if (tau_res_mode_str == "CONSTANT") {
        tau_res_mode = CONSTANT;
    } else if (tau_res_mode_str == "DISTRIBUTION") {
        tau_res_mode = DISTRIBUTION;
    } else {
        std::cerr <<
          "Invalid tau_res mode: " +
          tau_res_mode_str +
          ". Must be CONSTANT or DISTRIBUTION." << std::endl;
        throw(0);
    }
    std::cout << "> Conditions.tau_res_Mode = " << tauResModeString(tau_res_mode) << std::endl;
    switch(tau_res_mode) {
        case CONSTANT: {
            tau_res_constant = config->get_qualified_as<double>("Conditions.tau_res").value_or(DEFAULT_TAU_RES_CONSTANT);
            std::cout << "> Conditions.tau_res = " << tau_res_constant << std::endl;
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

    // Initialization
    // TODO - various initialization methods

    // Output
    check_interval = config->get_qualified_as<unsigned int>("Output.check_interval").value_or(DEFAULT_CHECK_INTERVAL);
    std::cout << "> Output.check_interval = " << check_interval << std::endl;
    std::vector<std::string> check_variable_names_ = config->get_qualified_array_of<std::string>("Output.check_variable_names").value_or(DEFAULT_CHECK_VARIABLE_NAMES);
    if (check_variable_names_.size() > 0) {
        check_variable_names.assign(check_variable_names_.begin(), check_variable_names_.end());
        std::cout << "> Output.check_variable_names: ";
        std::for_each(check_variable_names.begin(), check_variable_names.end(), [](std::string n){ std::cout << n << ", "; });
        std::cout << std::endl;
    }
    check_verbose = config->get_qualified_as<bool>("Output.check_verbose").value_or(DEFAULT_CHECK_VERBOSE);
    std::cout << "> Output.check_verbose = " << check_verbose << std::endl;
    write_interval = config->get_qualified_as<unsigned int>("Output.write_interval").value_or(DEFAULT_WRITE_INTERVAL);
    if (write_interval < 0) write_interval = check_interval;
    std::cout << "> Output.write_interval = " << write_interval << std::endl;

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
            case CONSTANT: {
                dt_step = 0.1 * std::min(tau_res_constant, tau_mix);
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
        std::cout << "Setting dt automatically: " << dt_step << std::endl;
    }
    if (dt_sub_target < 0.0) {
        dt_sub_target = 0.04 * tau_mix;
        std::cout << "Setting dt_sub automatically: " << dt_sub_target << std::endl;
    }

    // Create distributions and random engines for each thread
    for (int it = 0; it < omp_get_max_threads(); it++) {
        std::uniform_int_distribution<unsigned int> uni_int(0, n_particles);
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
    xtemp1.resize(n_state_variables);
    xtemp2.resize(n_state_variables);
    xmean_old.resize(n_state_variables);
    xvar_old.resize(n_state_variables);
    Y_fuel.resize(n_species);
    Y_ox.resize(n_species);
    Y_phi.resize(n_species);

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
    gasvec[0]->getMassFractions(Y_phi.data());
    double Zeq = gasvec[0]->mixtureFraction(comp_fuel, comp_ox);

    // // Compute effective C and H numbers in fuel
    // gasvec[0]->setState_TPX(T_fuel, P, comp_fuel);
    // double C_eff = 0.0;
    // double H_eff = 0.0;
    // int iC = gasvec[0]->elementIndex("C");
    // int iH = gasvec[0]->elementIndex("H");
    // for (size_t k = 0; k < n_species; k++) {
    //     C_eff += gasvec[0]->moleFraction(k) * gasvec[0]->nAtoms(k, iC);
    //     H_eff += gasvec[0]->moleFraction(k) * gasvec[0]->nAtoms(k, iH);
    // }

    // // Compute stoichiometric fuel/air ratio
    // double s = 32.0 * (C_eff + (H_eff/4.0)) / (12.0*C_eff + H_eff);
    // double YFF = 0.0;
    // for (int k = 0; k < n_species; k++) {
    //     if (Y_fuel[k] > 1.0e-3) {
    //         YFF += Y_fuel[k];
    //     }
    // }
    // double YOO = Y_ox[gasvec[0]->speciesIndex("O2")];
    // double S = s * YFF / YOO;

    // Compute equilibrium state (for initialization)
    gasvec[0]->setState_TPY(T_init, P, Y_phi.data());
    gasvec[0]->equilibrate("HP");
    double T_equil = gasvec[0]->temperature();
    double h_equil = gasvec[0]->enthalpy_mass();
    std::vector<double> Y_equil(n_species);
    gasvec[0]->getMassFractions(Y_equil.data());

    // Initialize particles
    std::cout << "Initializing particles..." << std::endl;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for
        for (int ip = 0; ip < n_particles * n_stat; ip++) {
            pvec[ip].setP(&P);
            pvec[ip].setIndex(ip);
            pvec[ip].setnSpecies(n_species);
            pvec[ip].setMass(1.0); // TODO: check this (fuel particles lighter?)
            pvec[ip].seth(h_equil);
            pvec[ip].setY(Y_equil.data());
            switch (tau_res_mode) {
                case CONSTANT: {
                    pvec[ip].setAge(dists_uni_real[tid](rand_engines[tid]) * tau_res_constant);
                    pvec[ip].setTauRes(tau_res_constant);
                    break;
                }
                case DISTRIBUTION: {
                    pvec[ip].setAge(0.0);
                    pvec[ip].setTauRes(tau_res_hist.rand(dists_uni_real[tid], rand_engines[tid]));
                    break;
                }
            }
            // if (ip <= (n_particles-1) / 2.0) {
            //     pvec[ip].setAge(0.0);
            //     pvec[ip].seth(h_fuel);
            //     pvec[ip].setY(Y_fuel.data());
            // } else {
            //     pvec[ip].setAge(1.0);
            //     pvec[ip].seth(h_ox);
            //     pvec[ip].setY(Y_ox.data());
            // }
            // pvec[ip].print();
        }
    }
    std::cout << "Done initializing particles." << std::endl;

    // Initialize injectors
    for (int iinj = 0; iinj < 2; iinj++) {
        injvec.push_back(Injector(iinj, n_species));
    }

    // Fuel injector
    injvec[0].seth(h_fuel);
    injvec[0].setY(Y_fuel.data());
    injvec[0].setFlow(Zeq);

    // Oxidizer injector
    injvec[1].seth(h_ox);
    injvec[1].setY(Y_ox.data());
    injvec[1].setFlow(1-Zeq);

    for (Injector& inj : injvec) {
        inj.print();
    }

    // Initialize statistics
    meanState(&xmean_old, true, true);
    varianceState(&xvar_old, true, true);

    // Initialize variable functions
    for (int iv = 0; iv < n_state_variables; iv++) {
        variable_functions.push_back(
            [this, iv](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
                return pvec[ip].state(iv); });
    }

    // Initialize derived variables
    derived_variable_names.push_back("age");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getAge(); });
    n_derived_variables++;

    derived_variable_names.push_back("tau_res");
    variable_functions.push_back(
        [this](std::shared_ptr<Cantera::ThermoPhase> gas, int ip) {
            return pvec[ip].getTauRes(); });
    n_derived_variables++;

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
}

void PartiallyStirredReactor::run() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Begin time stepping..." << std::endl;
    while (!runDone()) {
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

void PartiallyStirredReactor::check() {
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
    std::cout << "> rerror = " << rerror << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
}

std::string PartiallyStirredReactor::variableName(int iv) {
    if (iv == c_offset_h) {
        return "h";
    } else if (iv < n_state_variables) {
        std::string prefix = "Y_";
        return prefix + gasvec[0]->speciesName(iv-c_offset_Y);
    } else {
        return derived_variable_names[iv - n_state_variables];
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
    if (check_interval > 0) {
        if ((step % check_interval) == 0) {
            check();
        }
    }

    // Copy state to other stats region
    copyState();

    // Adjust time step
    if ((t_stop > 0) && ((t + dt_step) > t_stop)) {
        dt_step = t_stop - t;
    }

    int n_substeps = 1 + std::round(dt_step / dt_sub_target);
    double dt_sub = dt_step / n_substeps;

    // Take substeps
    // for (Particle& p : pvec) p.print();
    subStepInflow(dt_step);
    // for (Particle& p : pvec) p.print();
    for (int isub = 0; isub < n_substeps; isub++) {
        subStepMix(dt_sub);
        // for (Particle& p : pvec) p.print();
        subStepReact(dt_sub);
        // for (Particle& p : pvec) p.print();
    }
    // Increment counters
    incrementAge();
    step++;
    t += dt_step;

    // Calculate convergence metric
    calcConvergence();
}

void PartiallyStirredReactor::subStepInflow(double dt) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for
        // Iterate over particles and recycle
        for (int ip = 0; ip < n_particles; ip++) {
            if (pvec[ip].tooOld()) {
                double p_inj = dists_uni_real[tid](rand_engines[tid]); // Probability: which injector to use for new particle
                recycleParticle(ip, p_inj, tid);
            }
        }
    }
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
            throw Cantera::NotImplementedError("PartiallyStirredReactor::subStepMix",
                                               "CURL not implemented.");
            break;
        }
        case MOD_CURL: {
            // TODO - currently assumes equal weight particles
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
                    unsigned int ip1 = dists_uni_int[tid](rand_engines[tid]);
                    unsigned int ip2 = dists_uni_int[tid](rand_engines[tid]);
                    double a = dists_uni_real[tid](rand_engines[tid]);
                    pvec[ip1] += (1.0/2.0) * a * (pvec[ip2] - pvec[ip1]);
                    pvec[ip2] += (1.0/2.0) * a * (pvec[ip1] - pvec[ip2]);
                }
            }
            break;
        }
        case IEM: {
            meanState(&xtemp1, false, true);
            Particle pmean;
            pmean.setnSpecies(n_species);
            pmean.setState(xtemp1.data());
#pragma omp parallel for
            for (int ip = 0; ip < n_particles; ip++) {
                pvec[ip] += -(1.0/2.0) * (dt/tau_mix) * (pvec[ip] - pmean);
            }
            break;
        }
        case EMST: {
            throw Cantera::NotImplementedError("PartiallyStirredReactor::subStepMix",
                                               "EMST not implemented.");
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
    pvec[ip].setAge(0.0);
    switch (tau_res_mode) {
        case CONSTANT: {
            pvec[ip].setTauRes(tau_res_constant);
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
    // std::cout << "Recycling particle " << ip
    //     << " to injector " << iinj
    //     << " (p_inj = " << p_inj << ")" << std::endl;
    // pvec[ip].print();
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
    if (rerror <= rtol && (step >= min_steps_converge)) {
        std::cout << "Reached termination condition: rerror (" <<
            rerror << ") <= rtol (" << rtol << ") at step " << step << std::endl;
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

void PartiallyStirredReactor::writeStats() {

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
