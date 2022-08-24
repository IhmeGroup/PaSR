#include "Particle.h"

Particle::Particle() {

}

Particle::~Particle() {

}

void Particle::react(const double dt) {
    throw Cantera::NotImplementedError("Particle::react");
}