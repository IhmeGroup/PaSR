#include "Particle.h"

Particle::Particle() {

}

void Particle::print() {
    std::cout << "a: " << a() << "\t h: " << h() << std::endl;
}

void Particle::react(const double dt) {
    throw Cantera::NotImplementedError("Particle::react");
}

Particle::~Particle() {
}