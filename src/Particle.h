#include <iostream>

#include "cantera/base/ctexceptions.h"

const size_t c_offset_a = 0; // age
const size_t c_offset_T = 1; // temperature
const size_t c_offset_Y = 2; // mass fraction

class Particle {
public:
    explicit Particle();
    ~Particle();

    double& a() {
        return solvec[index(c_offset_a)];
    }

    double& T() {
        return solvec[index(c_offset_T)];
    }

    double& Y(const int& k) {
        return solvec[index(c_offset_Y + k)];
    }

    void setSolVec(double* solvec_) {
        solvec = solvec_;
    }

    void setOffset(const int& offset_) {
        offset = offset_;
    }

    void setnsp(const int& nsp_) {
        nsp = nsp_;
    }

    void setnv(const int& nv_) {
        nv = nv_;
    }

    void seta(const double& a) {
        solvec[index(c_offset_a)] = a;
    }

    void setT(const double& T) {
        solvec[index(c_offset_T)] = T;
    }

    void setY(const double* Y) {
        for (int k = 0; k < nsp; k++) {
            solvec[index(c_offset_Y + k)] = Y[k];
        }
    }

    void setState(const double& a, const double& T, double* Y) {
        seta(a);
        setT(T);
        setY(Y);
    }

    void setState(const double* state);

    void print();

    void react(const double dt);

protected:
    const int index(const int& i) {
        return (offset*nv) + i;
    }

    double* solvec;
    int offset;
    int nsp;
    int nv;

private:
};
