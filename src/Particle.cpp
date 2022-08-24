#include "Particle.h"

Particle::Particle() {

}

Particle::~Particle() {

}

void Particle::print() {
    std::cout << "a: " << a() << "\t T: " << T() << std::endl;
}

void Particle::react(const double dt) {
    throw Cantera::NotImplementedError("Particle::react");
}