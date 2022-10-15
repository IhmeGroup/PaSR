#include "Particle.h"

Particle::Particle() :
    id(0)
{
    setnSpecies(0);
}

Particle::Particle(const int& id, const unsigned int& n_species) {
    setID(id);
    setnSpecies(n_species);
}

void Particle::print(double threshold, std::shared_ptr<Cantera::ThermoPhase> gas) {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Particle " << id << ":" << std::endl;
    std::cout << "> mass: " << getMass() << std::endl;
    std::cout << "> age: " << getAge() << std::endl;
    if (gas) {
        std::cout << "> T: " << T(gas) << std::endl;
    }
    std::cout << "> h: " << h() << std::endl;
    std::cout << "> Y: " << std::endl;
    for (int k = 0; k < n_species; k++) {
        if (std::abs(Y(k)) >= threshold) {
            std::cout << ">   " << k << ": " << Y(k) << std::endl;
        }
    }
    std::cout << "--------------------------------------------------" << std::endl;
}

void Particle::react(Cantera::ReactorNet* rnet, double dt) {
    Cantera::Reactor* reactor = &rnet->reactor(0);
    Cantera::ThermoPhase* gas = &reactor->contents();
    gas->setState_PY(P(), Y());
    gas->setState_HP(h(), P());
    reactor->syncState();
    reactor->setInitialVolume(mass / gas->density());
    rnet->advance(rnet->time() + dt);
    seth(gas->enthalpy_mass());
    setY(gas->massFractions());
}

Particle::~Particle() {
}