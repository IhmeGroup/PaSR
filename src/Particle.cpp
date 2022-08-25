#include "Particle.h"

Particle::Particle() {

}

void Particle::initialize(const std::string& mech_filename) {
    sol = Cantera::newSolution(mech_filename);
    gas = sol->thermo();
    reactor = new Cantera::IdealGasConstPressureReactor();
    reactor->insert(sol);
    rnet = new Cantera::ReactorNet();
    rnet->addReactor(*reactor);
}

void Particle::print() {
    std::cout << "a: " << a() << "\t h: " << h() << std::endl;
}

void Particle::react(const std::string& mech_filename, const double& dt) {
    gas->setState_PY(P(), Y());
    gas->setState_HP(h(), P());
    reactor->setInitialVolume(mass / gas->density());
    rnet->advance(rnet->time() + dt);
    seth(gas->enthalpy_mass());
    setY(gas->massFractions());
}

Particle::~Particle() {
}