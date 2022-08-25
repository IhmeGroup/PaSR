#include "Injector.h"

Injector::Injector() :
    index(0)
{

}

Injector::Injector(const int& index, const unsigned int& nsp) {
    setIndex(index);
    setnsp(nsp);
}

void Injector::print() {
    std::cout << "----------" << std::endl;
    std::cout << "Injector " << index << ":" << std::endl;
    std::cout << "> h: " << h() << std::endl;
    std::cout << "> Y: " << std::endl;
    for (int k = 0; k < nsp; k++) {
        if (Y(k) > 1.0e-3) {
            std::cout << ">   " << k << ": " << Y(k) << std::endl;
        }
    }
    std::cout << "----------" << std::endl;
}

Injector::~Injector() {
    
}