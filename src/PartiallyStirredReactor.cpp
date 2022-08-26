#include <math.h>
#include <algorithm>
#include <omp.h>

#include "../cpptoml/include/cpptoml.h"

#include "PartiallyStirredReactor.h"

#pragma omp declare \
    reduction( \
        vec_double_plus : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            std::plus<double>())) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare \
    reduction( \
        vec_double_min : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            [](double x1, double x2) { return std::min(x1, x2); })) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare \
    reduction( \
        vec_double_max : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            [](double x1, double x2) { return std::max(x1, x2); })) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))

PartiallyStirredReactor::PartiallyStirredReactor(const std::string& input_filename_) :
    input_filename(input_filename_),
    step(0), t(0.0), p_out(0.0)
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
    np = config->get_qualified_as<unsigned int>("Numerics.n_particles").value_or(DEFAULT_N_PARTICLES);
    std::cout << "> Numerics.n_particles = " << np << std::endl;
    n_steps = config->get_qualified_as<unsigned int>("Numerics.n_steps").value_or(DEFAULT_N_STEPS);
    if (n_steps > 0)
        std::cout << "> Numerics.n_steps = " << n_steps << std::endl;
    t_stop = config->get_qualified_as<double>("Numerics.t_stop").value_or(DEFAULT_T_STOP);
    if (t_stop > 0)
        std::cout << "> Numerics.t_stop = " << t_stop << std::endl;
    dt_step = config->get_qualified_as<double>("Numerics.dt").value_or(DEFAULT_DT);
    if (dt_step > 0)
        std::cout << "> Numerics.dt = " << dt_step << std::endl;
    dt_sub_target = config->get_qualified_as<double>("Numerics.dt_sub").value_or(DEFAULT_DT_SUB);
    if (dt_sub_target > 0)
        std::cout << "> Numerics.dt_sub = " << dt_sub_target << std::endl;

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
    tau_res = config->get_qualified_as<double>("Conditions.tau_res").value_or(DEFAULT_TAU_RES);
    std::cout << "> Conditions.tau_res = " << tau_res << std::endl;
    tau_mix = config->get_qualified_as<double>("Conditions.tau_mix").value_or(DEFAULT_TAU_MIX);
    std::cout << "> Conditions.tau_mix = " << tau_mix << std::endl;

    // Output
    check_interval = config->get_qualified_as<unsigned int>("Output.check_interval").value_or(DEFAULT_CHECK_INTERVAL);
    std::cout << "> Output.check_interval = " << check_interval << std::endl;

    std::cout << "--------------------------------------------------" << std::endl;
}

void PartiallyStirredReactor::initialize() {
    step = 0;
    t = 0.0;
    p_out = 0.0;
    p_mix = 0.0;

    // Set step sizes
    if (dt_step <= 0.0) {
        dt_step = 0.1 * std::min(tau_res, tau_mix);
        std::cout << "Setting dt automatically: " << dt_step << std::endl;
    }
    if (dt_sub_target < 0.0) {
        dt_sub_target = 0.04 * tau_mix;
        std::cout << "Setting dt_sub automatically: " << dt_sub_target << std::endl;
    }

    // Generate random seed for each thread
    seedvec.resize(omp_get_max_threads());
    for (unsigned int& seed : seedvec) {
        seed = rand();
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
    
    nsp = gasvec[0]->nSpecies();
    nv = nsp + 2;
    pvec.resize(np);
    xvec.resize(nv*np);
    xtemp.resize(nv);
    Y_fuel.resize(nsp);
    Y_ox.resize(nsp);
    Y_phi.resize(nsp);

    // Get compositions
    gasvec[0]->setState_TPX(T_fuel, P, comp_fuel);
    gasvec[0]->getMassFractions(Y_fuel.data());
    h_fuel = gasvec[0]->enthalpy_mass();
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
    // for (size_t k = 0; k < nsp; k++) {
    //     C_eff += gasvec[0]->moleFraction(k) * gasvec[0]->nAtoms(k, iC);
    //     H_eff += gasvec[0]->moleFraction(k) * gasvec[0]->nAtoms(k, iH);
    // }

    // // Compute stoichiometric fuel/air ratio
    // double s = 32.0 * (C_eff + (H_eff/4.0)) / (12.0*C_eff + H_eff);
    // double YFF = 0.0;
    // for (int k = 0; k < nsp; k++) {
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
    std::vector<double> Y_equil(nsp);
    gasvec[0]->getMassFractions(Y_equil.data());

    // Initialize particles
    std::cout << "Initializing particles..." << std::endl;
#pragma omp parallel for
    for (int ip = 0; ip < np; ip++) {
        pvec[ip].setP(&P);
        pvec[ip].setIndex(ip);
        pvec[ip].setnsp(nsp);
        pvec[ip].setMass(1.0); // TODO: check this
        pvec[ip].setAge(0.0);
        pvec[ip].seth(h_equil);
        pvec[ip].setY(Y_equil.data());
        // if (ip <= (np-1) / 2.0) {
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
    std::cout << "Done initializing particles." << std::endl;

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
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Begin time stepping..." << std::endl;
    while (!runDone()) {
        takeStep();
    }
    std::cout << "Done." << std::endl;
}

void PartiallyStirredReactor::print() {
    
}

void PartiallyStirredReactor::check() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Starting step: " << step << std::endl;
    std::cout << "> t: " << t << std::endl;
    std::cout << "> a: min = " << minAge() <<
        ", fmean = " << meanAge(true) <<
        ", max = " << maxAge() << std::endl;
}

void PartiallyStirredReactor::takeStep() {
    // Print status
    if (check_interval > 0) {
        if ((step % check_interval) == 0) {
            check();
        }
    }

    // Adjust time step
    if ((t_stop > 0) && ((t + dt_step) > t_stop)) {
        dt_step = t_stop - t;
    }

    int n_substeps = 1 + round(dt_step / dt_sub_target);
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
}

void PartiallyStirredReactor::subStepInflow(double dt) {
    p_out += np * dt / tau_res; // Fractional particle count to recycle
    int np_out = round(p_out); // Round to integer
    p_out -= np_out; // Hold on to remainder for next step

    unsigned int* s;
#pragma omp parallel private(s)
    {
        // Grab this thread's seed
        s = &seedvec[omp_get_thread_num()];
#pragma omp for
        // Iterate over particles and recycle
        for (int ip_out = 0; ip_out < np_out; ip_out++) {
            unsigned int ip = rand_r(s) % np; // Index of particle to recycle
            double p_inj = (rand_r(s) / (double)RAND_MAX); // Probability: which injector to use for new particle
            recycleParticle(ip, p_inj);
        }
        s++;
    }
}

void PartiallyStirredReactor::subStepMix(double dt) {
    switch(mixing_model) {
        case NO_MIX: {
            // Do nothing
            break;
        }
        case FULL_MIX: {
            meanState(&xtemp, true);
            // Set all particles to Favre mean state
#pragma omp parallel for
            for (int ip = 0; ip < np; ip++) {
                pvec[ip].setState(xtemp.data());
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
            p_mix += np * dt / tau_mix;
            int np_mix = round(p_mix);
            p_mix -= np_mix;
            unsigned int* s;
#pragma omp parallel private(s)
            {
                // Grab this thread's seed
                s = &seedvec[omp_get_thread_num()];
#pragma omp for
                // Iterate over pairs and mix
                for (int ipair = 0; ipair < np_mix; ipair++) {
                    unsigned int ip1 = rand_r(s) % np;
                    unsigned int ip2 = rand_r(s) % np;
                    double a = (rand_r(s) / (double)RAND_MAX);
                    pvec[ip1] += (1.0/2.0) * a * (pvec[ip2] - pvec[ip1]);
                    pvec[ip2] += (1.0/2.0) * a * (pvec[ip1] - pvec[ip2]);
                }
                s++;
            }
            break;
        }
        case IEM: {
            meanState(&xtemp, true);
            Particle pmean;
            pmean.setnsp(nsp);
            pmean.setState(xtemp.data());
#pragma omp parallel for
            for (int ip = 0; ip < np; ip++) {
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
            throw Cantera::NotImplementedError("PartiallyStirredReactor::subStepMix",
                                               "Invalid mixing model.");
        }
    }
}

void PartiallyStirredReactor::subStepReact(double dt) {
#pragma omp parallel for
    for (int ip = 0; ip < np; ip++) {
        // Use the current thread's reactor for calculation
        pvec[ip].react(rnetvec[omp_get_thread_num()], dt);
    }
}

void PartiallyStirredReactor::incrementAge() {
#pragma omp parallel for
    for (int ip = 0; ip < np; ip++) {
        pvec[ip].getAge() += dt_step;
    }
}

void PartiallyStirredReactor::recycleParticle(unsigned int ip, double p_inj) {
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
    pvec[ip].seth(injvec[iinj].h());
    pvec[ip].setY(injvec[iinj].Y().data());
    // std::cout << "Recycling particle " << ip
    //     << " to injector " << iinj
    //     << " (p_inj = " << p_inj << ")" << std::endl;
    // pvec[ip].print();
}

double PartiallyStirredReactor::mean(std::function<double(int)> xfunc, bool favre) {
    double rhosum = 0.0;
    double xsum = 0.0;
#pragma omp parallel for reduction(+:rhosum,xsum)
    for (int ip = 0; ip < np; ip++) {
        double rho = 0.0;
        if (favre) {
            rho = pvec[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        xsum += rho * xfunc(ip);
    }
    return xsum / rhosum;
}

double PartiallyStirredReactor::meanState(int iv, bool favre) {
    double rhosum = 0.0;
    double xsum = 0.0;
#pragma omp parallel for reduction(+:rhosum,xsum)
    for (int ip = 0; ip < np; ip++) {
        double rho = 0.0;
        if (favre) {
            rho = pvec[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        xsum += rho * pvec[ip].state(iv);
    }
    return xsum / rhosum;
}

void PartiallyStirredReactor::meanState(std::vector<double>* xmeanvec, bool favre) {

    // Sum across particles
    std::vector<double> xsumvec(xmeanvec->size(), 0.0);
    double rhosum = 0.0;
#pragma omp parallel for reduction(+:rhosum) reduction(vec_double_plus:xsumvec)
    for (int ip = 0; ip < np; ip++) {
        double rho;
        if (favre) {
            rho = pvec[ip].rho(gasvec[omp_get_thread_num()]);
        } else {
            rho = 1.0;
        }
        rhosum += rho;
        for (int iv = 0; iv < nv; iv++) {
            xsumvec[iv] += rho * pvec[ip].state(iv);
        }
    }

    // Divide by particle count and write to mean vec
    for (int iv = 0; iv < nv; iv++) {
        (*xmeanvec)[iv] = xsumvec[iv] / rhosum;
    }
}

double PartiallyStirredReactor::min(std::function<double(int)> xfunc) {
    double minval = xfunc(0);
#pragma omp parallel for reduction(min:minval)
    for (int ip = 0; ip < np; ip++) {
        minval = std::min(minval, xfunc(ip));
    }
    return minval;
}

double PartiallyStirredReactor::minState(int iv) {
    double minval = pvec[0].state(iv);
#pragma omp parallel for reduction(min:minval)
    for (int ip = 0; ip < np; ip++) {
        minval = std::min(minval, pvec[ip].state(iv));
    }
    return minval;
}

void PartiallyStirredReactor::minState(std::vector<double>* minvec) {

    // Initialize min with first particle
    std::vector<double> minvec_temp(minvec->size());
    for (int iv = 0; iv < nv; iv++) {
        minvec_temp[iv] = pvec[0].state(iv);
    }

    // Min across particles
#pragma omp parallel for reduction(vec_double_min:minvec_temp)
    for (int ip = 0; ip < np; ip++) {
        for (int iv = 0; iv < nv; iv++) {
            minvec_temp[iv] = std::min(minvec_temp[iv], pvec[ip].state(iv));
        }
    }

    // Write to minvec
    for (int iv = 0; iv < nv; iv++) {
        (*minvec)[iv] = minvec_temp[iv];
    }
}

double PartiallyStirredReactor::max(std::function<double(int)> xfunc) {
    double maxval = xfunc(0);
#pragma omp parallel for reduction(max:maxval)
    for (int ip = 0; ip < np; ip++) {
        maxval = std::max(maxval, xfunc(ip));
    }
    return maxval;
}

double PartiallyStirredReactor::maxState(int iv) {
    double maxval = pvec[0].state(iv);
#pragma omp parallel for reduction(max:maxval)
    for (int ip = 0; ip < np; ip++) {
        maxval = std::max(maxval, pvec[ip].state(iv));
    }
    return maxval;
}

void PartiallyStirredReactor::maxState(std::vector<double>* maxvec) {

    // Initialize max with first particle
    std::vector<double> maxvec_temp(maxvec->size());
    for (int iv = 0; iv < nv; iv++) {
        maxvec_temp[iv] = pvec[0].state(iv);
    }

    // Max across particles
#pragma omp parallel for reduction(vec_double_max:maxvec_temp)
    for (int ip = 0; ip < np; ip++) {
        for (int iv = 0; iv < nv; iv++) {
            maxvec_temp[iv] = std::max(maxvec_temp[iv], pvec[ip].state(iv));
        }
    }

    // Write to maxvec
    for (int iv = 0; iv < nv; iv++) {
        (*maxvec)[iv] = maxvec_temp[iv];
    }
}

double PartiallyStirredReactor::meanAge(bool favre) {
    std::function<double(int)> age = [this](int ip) { return pvec[ip].getAge(); };
    return mean(age, favre);
}

double PartiallyStirredReactor::minAge() {
    std::function<double(int)> age = [this](int ip) { return pvec[ip].getAge(); };
    return min(age);
}

double PartiallyStirredReactor::maxAge() {
    std::function<double(int)> age = [this](int ip) { return pvec[ip].getAge(); };
    return max(age);
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
