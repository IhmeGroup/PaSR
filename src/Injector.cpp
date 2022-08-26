#include "Injector.h"

Injector::Injector() :
    index(0)
{

}

Injector::Injector(const int& index, const unsigned int& nsp) {
    setIndex(index);
    setnsp(nsp);
}

void Injector::print(double threshold) {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Injector " << index << ":" << std::endl;
    std::cout << "> flow: " << getFlow() << std::endl;
    std::cout << "> h: " << h() << std::endl;
    std::cout << "> Y: " << std::endl;
    for (int k = 0; k < nsp; k++) {
        if (Y(k) >= threshold) {
            std::cout << ">   " << k << ": " << Y(k) << std::endl;
        }
    }
    std::cout << "--------------------------------------------------" << std::endl;
}

Injector::~Injector() {
    
}