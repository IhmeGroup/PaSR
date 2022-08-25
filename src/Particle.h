#include <iostream>

#include "cantera/base/ctexceptions.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/zerodim.h"

const size_t c_offset_a = 0; // age
const size_t c_offset_h = 1; // temperature
const size_t c_offset_Y = 2; // mass fraction

class Particle {
public:
    explicit Particle();
    ~Particle();

    // Particle operator+= (const Particle& rhs)& {
    //     for (int iv = 0; iv < nv; iv++) {
    //         state(iv) += rhs.state(iv);
    //     }
    //     return *this;
    // }

    // friend Particle operator+ (Particle lhs, const Particle& rhs) {
    //     lhs += rhs;
    //     return lhs;
    // }

    double& P() {
        return *m_P;
    }

    double* state() {
        return &xvec[solIndex(0)];
    }

    double& state(const int& i) {
        return xvec[solIndex(i)];
    }

    double& a() {
        return xvec[solIndex(c_offset_a)];
    }

    double& h() {
        return xvec[solIndex(c_offset_h)];
    }

    double& T() {
        throw Cantera::NotImplementedError("Particle::T");
    }

    double* Y() {
        return &xvec[solIndex(c_offset_Y)];
    }

    double& Y(const int& k) {
        return xvec[solIndex(c_offset_Y + k)];
    }

    void setSolVec(double* xvec_) {
        xvec = xvec_;
    }

    void setIndex(const int& index_) {
        index = index_;
    }

    void setnsp(const int& nsp_) {
        nsp = nsp_;
    }

    void setnv(const int& nv_) {
        nv = nv_;
    }
    
    void setMass(const double& mass_) {
        mass = mass_;
    }

    void setP(double* P_) {
        m_P = P_;
    }

    void seta(const double& a) {
        xvec[solIndex(c_offset_a)] = a;
    }

    void seth(const double& h) {
        xvec[solIndex(c_offset_h)] = h;
    }

    void setT(const double& T) {
        throw Cantera::NotImplementedError("Particle::setT");
    }

    void setY(const double* Y) {
        for (int k = 0; k < nsp; k++) {
            xvec[solIndex(c_offset_Y + k)] = Y[k];
        }
    }

    void setState(const double& a, const double& h, double* Y) {
        seta(a);
        seth(h);
        setY(Y);
    }

    void setState(const double* state) {
        throw Cantera::NotImplementedError("Particle::setState(const double* state)");
    }

    void print();

    void react(Cantera::ReactorNet* rnet, const double& dt);

protected:
    int solIndex(const int& i) {
        return (index*nv) + i;
    }

    double* xvec = nullptr;
    double* m_P = nullptr;
    int index;
    int nsp;
    int nv;
    double mass;

private:
};
