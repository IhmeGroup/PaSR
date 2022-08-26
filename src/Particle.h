#include <iostream>

#include "cantera/base/ctexceptions.h"
#include "cantera/core.h"
#include "cantera/zerodim.h"

const size_t c_offset_h = 0; // enthalpy
const size_t c_offset_Y = 1; // mass fraction

class Particle {
public:
    explicit Particle();
    ~Particle();

    Particle operator+= (const Particle& rhs)& {
        for (int iv = 0; iv < nv; iv++) {
            xvec[iv] += rhs.xvec[iv];
        }
        return *this;
    }

    friend Particle operator+ (Particle lhs, const Particle& rhs) {
        lhs += rhs;
        return lhs;
    }

    Particle operator-= (const Particle& rhs)& {
        for (int iv = 0; iv < nv; iv++) {
            xvec[iv] -= rhs.xvec[iv];
        }
        return *this;
    }

    friend Particle operator- (Particle lhs, const Particle& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Particle operator*= (const double& rhs)& {
        for (int iv = 0; iv < nv; iv++) {
            xvec[iv] *= rhs;
        }
        return *this;
    }

    friend Particle operator* (Particle lhs, const double& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Particle operator* (const double& lhs, Particle rhs) {
        rhs *= lhs;
        return rhs;
    }

    double& P() {
        return *m_P;
    }

    double* state() {
        return &xvec[0];
    }

    double& state(const int& i) {
        return xvec[i];
    }

    double& getAge() {
        return age;
    }

    double& h() {
        return xvec[c_offset_h];
    }

    double& T() {
        throw Cantera::NotImplementedError("Particle::T");
    }

    double* Y() {
        return &xvec[c_offset_Y];
    }

    double& Y(const int& k) {
        return xvec[c_offset_Y + k];
    }

    void setIndex(const int& index_) {
        index = index_;
    }

    void setnsp(const int& nsp_) {
        nsp = nsp_;
        nv = nsp + c_offset_Y;
        xvec.resize(nv);
    }
    
    void setMass(const double& mass_) {
        mass = mass_;
    }

    void setP(double* P_) {
        m_P = P_;
    }

    void setAge(const double& age_) {
        age = age_;
    }

    void seth(const double& h) {
        xvec[c_offset_h] = h;
    }

    void setT(const double& T) {
        throw Cantera::NotImplementedError("Particle::setT");
    }

    void setY(const double* Y) {
        for (int k = 0; k < nsp; k++) {
            xvec[c_offset_Y + k] = Y[k];
        }
    }

    void setState(const double& h, double* Y) {
        seth(h);
        setY(Y);
    }

    void setState(const double* state) {
        throw Cantera::NotImplementedError("Particle::setState(const double* state)");
    }

    void print();

    void react(Cantera::ReactorNet* rnet, const double& dt);

protected:
    int index;
    int nsp;
    int nv;
    double* m_P = nullptr;
    std::vector<double> xvec;
    double age;
    double mass;

private:
};
