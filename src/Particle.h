#include <iostream>

#include "cantera/base/ctexceptions.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/ReactorNet.h"

const int c_offset_h = 0; // enthalpy
const int c_offset_Y = 1; // mass fraction

class Particle {
public:
    explicit Particle();
    ~Particle();

    // Particle operator= (const Particle& rhs) {
    //     for (int iv = 0; iv < nv; iv++) {
    //         xvec[iv] = rhs.xvec[iv];
    //     }
    //     return *this;
    // }

    // Particle copyAll(const Particle& rhs) {
    //     for (int iv = 0; iv < nv; iv++) {
    //         xvec[iv] = rhs.xvec[iv];
    //     }
    //     index = rhs.index;
    //     m_P = rhs.m_P;
    //     age = rhs.age;
    //     mass = rhs.mass;
    //     tau_res = rhs.tau_res;
    //     n_recycles = rhs.n_recycles;
    //     return *this;
    // }

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

    double& state(int i) {
        return xvec[i];
    }

    double& getAge() {
        return age;
    }

    double& getTauRes() {
        return tau_res;
    }

    int& getnRecycles() {
        return n_recycles;
    }

    bool tooOld() {
        return age >= tau_res;
    }

    double& h() {
        return xvec[c_offset_h];
    }

    double* Y() {
        return &xvec[c_offset_Y];
    }

    double& Y(int k) {
        return xvec[c_offset_Y + k];
    }

    double rho(std::shared_ptr<Cantera::ThermoPhase> gas) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->density();
    }

    double rho(Cantera::ThermoPhase* gas) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->density();
    }

    double T(std::shared_ptr<Cantera::ThermoPhase> gas) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->temperature();
    }

    double T(Cantera::ThermoPhase* gas) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->temperature();
    }

    double Z(std::shared_ptr<Cantera::ThermoPhase> gas,
             const std::string& comp_fuel,
             const std::string& comp_ox) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->mixtureFraction(comp_fuel, comp_ox);
    }

    double Z(Cantera::ThermoPhase* gas,
             const std::string& comp_fuel,
             const std::string& comp_ox) {
        gas->setState_PY(P(), Y());
        gas->setState_HP(h(), P());
        return gas->mixtureFraction(comp_fuel, comp_ox);
    }

    void setIndex(int index_) {
        index = index_;
    }

    void setnsp(int nsp_) {
        nsp = nsp_;
        nv = nsp + c_offset_Y;
        xvec.resize(nv);
    }
    
    void setMass(double mass_) {
        mass = mass_;
    }

    void setP(double* P_) {
        m_P = P_;
    }

    void setAge(double age_) {
        age = age_;
    }
    
    void setTauRes(double tau_res_) {
        tau_res = tau_res_;
    }

    void seth(double h) {
        xvec[c_offset_h] = h;
    }

    void setY(const double* Y) {
        for (int k = 0; k < nsp; k++) {
            xvec[c_offset_Y + k] = Y[k];
        }
    }

    void setState(double h, const double* Y) {
        seth(h);
        setY(Y);
    }

    void setState(const double* state) {
        for (int iv = 0; iv < nv; iv++) {
            xvec[iv] = state[iv];
        }
    }

    void print(double threshold = 1.0e-14);

    void react(Cantera::ReactorNet* rnet, double dt);

protected:
    int index;
    int nsp;
    int nv;
    double* m_P = nullptr;
    std::vector<double> xvec;
    double age;
    double mass;
    double tau_res;
    int n_recycles;

private:
};
