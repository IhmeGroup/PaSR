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

    void initialize(const std::string& mech_filename);

    double& P() {
        return *m_P;
    }

    double& a() {
        return solvec[solIndex(c_offset_a)];
    }

    double& h() {
        return solvec[solIndex(c_offset_h)];
    }

    double& T() {
        throw Cantera::NotImplementedError("Particle::T");
    }

    double* Y() {
        return &solvec[solIndex(c_offset_Y)];
    }

    double& Y(const int& k) {
        return solvec[solIndex(c_offset_Y + k)];
    }

    void setSolVec(double* solvec_) {
        solvec = solvec_;
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
        solvec[solIndex(c_offset_a)] = a;
    }

    void seth(const double& h) {
        solvec[solIndex(c_offset_h)] = h;
    }

    void setT(const double& T) {
        throw Cantera::NotImplementedError("Particle::setT");
    }

    void setY(const double* Y) {
        for (int k = 0; k < nsp; k++) {
            solvec[solIndex(c_offset_Y + k)] = Y[k];
        }
    }

    void setState(const double& a, const double& h, double* Y) {
        seta(a);
        seth(h);
        setY(Y);
    }

    void setState(const double* state);

    void print();

    void react(const std::string& mech_filename, const double& dt);

protected:
    int solIndex(const int& i) {
        return (index*nv) + i;
    }

    double* solvec = nullptr;
    double* m_P = nullptr;
    std::shared_ptr<Cantera::Solution> sol = nullptr;
    std::shared_ptr<Cantera::ThermoPhase> gas = nullptr;
    Cantera::IdealGasConstPressureReactor* reactor;
    Cantera::ReactorNet* rnet;
    int index;
    int nsp;
    int nv;
    double mass;

private:
};
