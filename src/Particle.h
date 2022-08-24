#include <iostream>
#include <vector>

#include "cantera/base/ctexceptions.h"

class Particle {
public:
    explicit Particle();
    ~Particle();

    double a();
    double T();
    double Y(int k);

    void setState(const double& T, const std::vector<double>* Y);
    void setState(const std::vector<double>* state);

    void react(const double dt);

protected:
    std::vector<double>* data;
    int offset;

private:
};
