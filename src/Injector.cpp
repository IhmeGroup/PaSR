#include "Injector.h"

Injector::Injector() :
    id(0)
{
    setnSpecies(0);
}

Injector::Injector(const int& id, const unsigned int& n_species) {
    setID(id);
    setnSpecies(n_species);
}

void Injector::print(double threshold) {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Injector " << id << ":" << std::endl;
    std::cout << "> flow: " << getFlow() << std::endl;
    std::cout << "> h: " << h() << std::endl;
    std::cout << "> Y: " << std::endl;
    for (int k = 0; k < n_species; k++) {
        if (Y(k) >= threshold) {
            std::cout << ">   " << k << ": " << Y(k) << std::endl;
        }
    }
    std::cout << "--------------------------------------------------" << std::endl;
}

Injector::~Injector() {
    
}