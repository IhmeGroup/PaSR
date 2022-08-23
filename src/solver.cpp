#include <omp.h>
#include "../external/tomlplusplus/toml.hpp"

#include "solver.h"

Solver::Solver(const std::string& input_filename_) :
    input_filename(input_filename_),
    step(0), t(0.0)
{
    parseInput();
}

void Solver::parseInput() {
    auto config = toml::parse_file(input_filename);

    // Mechanism
    mech_filename = config["Mechanism"]["name"].value_or("");
    std::cout << "Mechanism.name = " << mech_filename << std::endl;

    // Numerics
    n_particles = config["Numerics"]["n_particles"].value_or(0);
    std::cout << "Numerics.n_particles = " << n_particles << std::endl;
    n_steps = config["Numerics"]["n_steps"].value_or(-1);
    std::cout << "Numerics.n_steps = " << n_steps << std::endl;
    t_stop = config["Numerics"]["t_stop"].value_or(-1.0);
    std::cout << "Numerics.t_stop = " << t_stop << std::endl;
    dt = config["Numerics"]["dt"].value_or(0.0);
    std::cout << "Numerics.dt = " << dt << std::endl;

    // Models
    std::string mixing_model_str = config["Models"]["mixing_model"].value_or("FULL_MIX");
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
    std::cout << "Models.mixing_model = " << mixingModelString(mixing_model) << std::endl;

    // Conditions
    P = config["Conditions"]["pressure"].value_or(101325.0);
    std::cout << "Conditions.pressure = " << P << std::endl;
    comp_fuel = config["Conditions"]["comp_fuel"].value_or("");
    std::cout << "Conditions.comp_fuel = " << comp_fuel << std::endl;
    comp_ox = config["Conditions"]["comp_ox"].value_or("");
    std::cout << "Conditions.comp_ox = " << comp_ox << std::endl;
    T_fuel = config["Conditions"]["T_fuel"].value_or(300.0);
    std::cout << "Conditions.T_fuel = " << T_fuel << std::endl;
    T_ox = config["Conditions"]["T_ox"].value_or(300.0);
    std::cout << "Conditions.T_ox = " << T_ox << std::endl;
    phi_global = config["Conditions"]["phi_global"].value_or(1.0);
    std::cout << "Conditions.phi_global = " << phi_global << std::endl;
    tau_res = config["Conditions"]["tau_res"].value_or(1.0);
    std::cout << "Conditions.tau_res = " << tau_res << std::endl;
    tau_mix = config["Conditions"]["tau_mix"].value_or(1.0);
    std::cout << "Conditions.tau_mix = " << tau_mix << std::endl;

    // Computation
    n_threads = config["Computation"]["n_threads"].value_or(omp_get_max_threads());
    std::cout << "Computation.n_threads = " << n_threads << std::endl;
}

void Solver::run() {
    while (!runDone()) {
        takeStep();
    }
    std::cout << "Done." << std::endl;
}

void Solver::print() {
    std::cout << "step: " << step << "\tt: " << t << std::endl;
}

void Solver::takeStep() {
    print();
    step++;
    t += dt;
}

bool Solver::runDone() {
    if ((n_steps > 0) && (step >= n_steps)) {
        std::cout << "Reached termination condition: step >= n_steps" << std::endl;
        return true;
    }
    if ((t_stop > 0.0) && (t >= t_stop)) {
        std::cout << "Reached termination condition: t >= t_stop" << std::endl;
        return true;
    }
    return false;
}