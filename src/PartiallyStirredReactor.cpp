#include <math.h>
#include <omp.h>

#include "../cpptoml/include/cpptoml.h"

#include "PartiallyStirredReactor.h"

PartiallyStirredReactor::PartiallyStirredReactor(const std::string& input_filename_) :
    input_filename(input_filename_),
    step(0), t(0.0), p_out(0.0)
{
    parseInput();
}

void PartiallyStirredReactor::parseInput() {
    auto config = cpptoml::parse_file(input_filename);

    std::cout << "----------" << std::endl;
    std::cout << "Input file parameters:" << std::endl;

    // Mechanism
    mech_filename = config->get_qualified_as<std::string>("Mechanism.name").value_or("");
    std::cout << "> Mechanism.name = " << mech_filename << std::endl;

    // Numerics
    np = config->get_qualified_as<double>("Numerics.n_particles").value_or(0);
    std::cout << "> Numerics.n_particles = " << np << std::endl;
    n_steps = config->get_qualified_as<int>("Numerics.n_steps").value_or(-1);
    if (n_steps > 0)
        std::cout << "> Numerics.n_steps = " << n_steps << std::endl;
    t_stop = config->get_qualified_as<double>("Numerics.t_stop").value_or(-1.0);
    if (t_stop > 0)
        std::cout << "> Numerics.t_stop = " << t_stop << std::endl;
    dt = config->get_qualified_as<double>("Numerics.dt").value_or(-1.0);
    if (dt > 0)
        std::cout << "> Numerics.dt = " << dt << std::endl;

    // Models
    std::string mixing_model_str = config->get_qualified_as<std::string>("Models.mixing_model").value_or("FULL_MIX");
    if (mixing_model_str == "NO_MIX") {
        mixing_model = NO_MIX;
    }else if (mixing_model_str == "FULL_MIX") {
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
    P = config->get_qualified_as<double>("Conditions.pressure").value_or(101325.0);
    std::cout << "> Conditions.pressure = " << P << std::endl;
    comp_fuel = config->get_qualified_as<std::string>("Conditions.comp_fuel").value_or("");
    std::cout << "> Conditions.comp_fuel = " << comp_fuel << std::endl;
    comp_ox = config->get_qualified_as<std::string>("Conditions.comp_ox").value_or("");
    std::cout << "> Conditions.comp_ox = " << comp_ox << std::endl;
    T_fuel = config->get_qualified_as<double>("Conditions.T_fuel").value_or(300.0);
    std::cout << "> Conditions.T_fuel = " << T_fuel << std::endl;
    T_ox = config->get_qualified_as<double>("Conditions.T_ox").value_or(300.0);
    std::cout << "> Conditions.T_ox = " << T_ox << std::endl;
    T_init = config->get_qualified_as<double>("Conditions.T_init").value_or(1500.0);
    std::cout << "> Conditions.T_init = " << T_init << std::endl;
    phi_global = config->get_qualified_as<double>("Conditions.phi_global").value_or(1.0);
    std::cout << "> Conditions.phi_global = " << phi_global << std::endl;
    tau_res = config->get_qualified_as<double>("Conditions.tau_res").value_or(1.0);
    std::cout << "> Conditions.tau_res = " << tau_res << std::endl;
    tau_mix = config->get_qualified_as<double>("Conditions.tau_mix").value_or(1.0);
    std::cout << "> Conditions.tau_mix = " << tau_mix << std::endl;

    std::cout << "----------" << std::endl;
}

void PartiallyStirredReactor::initialize() {
    step = 0;
    t = 0.0;
    p_out = 0.0;

    if (dt <= 0.0) {
        dt = 0.1 * std::min(tau_res, tau_mix);
        std::cout << "Setting dt automatically: " << dt << std::endl;
    }

    sol = Cantera::newSolution(mech_filename);
    gas = sol->thermo();
    
    nsp = gas->nSpecies();
    nv = nsp + 2;
    pvec.resize(np);
    solvec.resize(nv*np);
    Y_fuel.resize(nsp);
    Y_ox.resize(nsp);
    Y_phi.resize(nsp);

    // Get compositions
    gas->setState_TPX(T_fuel, P, comp_fuel);
    gas->getMassFractions(Y_fuel.data());
    h_fuel = gas->enthalpy_mass();
    gas->setState_TPX(T_ox, P, comp_ox);
    gas->getMassFractions(Y_ox.data());
    h_ox = gas->enthalpy_mass();
    gas->setEquivalenceRatio(phi_global, comp_fuel, comp_ox);
    gas->getMassFractions(Y_phi.data());
    double Zeq = gas->mixtureFraction(comp_fuel, comp_ox);

    // // Compute effective C and H numbers in fuel
    // gas->setState_TPX(T_fuel, P, comp_fuel);
    // double C_eff = 0.0;
    // double H_eff = 0.0;
    // int iC = gas->elementIndex("C");
    // int iH = gas->elementIndex("H");
    // for (size_t k = 0; k < nsp; k++) {
    //     C_eff += gas->moleFraction(k) * gas->nAtoms(k, iC);
    //     H_eff += gas->moleFraction(k) * gas->nAtoms(k, iH);
    // }

    // // Compute stoichiometric fuel/air ratio
    // double s = 32.0 * (C_eff + (H_eff/4.0)) / (12.0*C_eff + H_eff);
    // double YFF = 0.0;
    // for (int k = 0; k < nsp; k++) {
    //     if (Y_fuel[k] > 1.0e-3) {
    //         YFF += Y_fuel[k];
    //     }
    // }
    // double YOO = Y_ox[gas->speciesIndex("O2")];
    // double S = s * YFF / YOO;

    // Compute equilibrium state (for initialization)
    gas->setState_TPY(T_init, P, Y_phi.data());
    gas->equilibrate("HP");
    double T_equil = gas->temperature();
    double h_equil = gas->enthalpy_mass();
    std::vector<double> Y_equil(nsp);
    gas->getMassFractions(Y_equil.data());

    // Initialize particles
    std::cout << "Initializing particles" << std::endl;
    for (int ip = 0; ip < np; ip++) {
        pvec[ip].setSolVec(solvec.data());
        pvec[ip].setP(&P);
        pvec[ip].setIndex(ip);
        pvec[ip].setnsp(nsp);
        pvec[ip].setnv(nv);
        pvec[ip].setMass(1.0); // TODO: check this
        pvec[ip].seta(0.0);
        pvec[ip].seth(h_equil);
        pvec[ip].setY(Y_equil.data());
        pvec[ip].initialize(mech_filename);
        // pvec[ip].print();
    }
    std::cout << "Done initializing particles" << std::endl;

    // Initialize injectors
    for (int iinj = 0; iinj < 2; iinj++) {
        injvec.push_back(Injector(iinj, nsp));
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
}

void PartiallyStirredReactor::run() {
    std::cout << "----------" << std::endl;
    std::cout << "Begin time stepping" << std::endl;
    while (!runDone()) {
        takeStep();
    }
    std::cout << "Done." << std::endl;
}

void PartiallyStirredReactor::print() {
    std::cout << "step: " << step << "\tt: " << t << std::endl;
}

void PartiallyStirredReactor::takeStep() {

    // Adjust time step
    if ((t_stop > 0) && ((t + dt) > t_stop)) {
        dt = t_stop - t;
    }

    // ====================
    // Inflow substep
    // ====================
    p_out += np * dt / tau_res; // Fractional particle count to recycle
    int np_out = round(p_out); // Round to integer
    p_out -= np_out; // Hold on to remainder for next step

    // Generate random seed for each thread
    std::vector<int> seedvec(omp_get_max_threads());
    for (int& seed : seedvec) {
        seed = rand();
    }

    unsigned int s;
#pragma omp parallel private(s)
    {
        // Grab this thread's seed
        s = seedvec[omp_get_thread_num()];
#pragma omp for
        // Iterate over particles and recycle
        for (int ip_out = 0; ip_out < np_out; ip_out++) {
            unsigned int ip = rand_r(&s) % np; // Index of particle to recycle
            double p_inj = (rand_r(&s) / (double)RAND_MAX); // Probability: which injector to use for new particle
            // std::cout << "thread: " << omp_get_thread_num() << " ip: " << ip << " p_inj: " << p_inj << std::endl;
            recycleParticle(ip, p_inj);
        }
    }

    // ====================
    // Mixing substep
    // ====================

    // ====================
    // Reaction substep
    // ====================
#pragma omp parallel for
    for (int ip = 0; ip < np; ip++) {
        pvec[ip].react(mech_filename, dt);
    }

    // Print status
    print();

    // Increment
    step++;
    t += dt;
}

void PartiallyStirredReactor::recycleParticle(const unsigned int& ip, const double& p_inj) {
    int iinj;
    double flow_sum = 0.0;
    for (int i = 0; i < injvec.size(); i++) {
        flow_sum += injvec[i].getFlow();
        if (p_inj <= flow_sum) {
            iinj = i;
            break;
        }
    }
    pvec[ip].seta(0.0);
    pvec[ip].seth(injvec[iinj].h());
    pvec[ip].setY(injvec[iinj].Y().data());
    // std::cout << "Recycling particle " << ip
    //     << " to injector " << iinj
    //     << " (p_inj = " << p_inj << ")" << std::endl;
    // pvec[ip].print();
}

bool PartiallyStirredReactor::runDone() {
    if ((n_steps > 0) && (step >= n_steps)) {
        std::cout << "Reached termination condition: step >= n_steps" << std::endl;
        return true;
    }
    if ((t_stop > 0.0) && (t >= t_stop)) {
        std::cout << "Reached termination condition: t >= t_stop" << std::endl;
        return true;
    }
    // TODO - residual-based termination
    return false;
}

PartiallyStirredReactor::~PartiallyStirredReactor() {
}
