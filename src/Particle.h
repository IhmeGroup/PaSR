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
    explicit Particle(const int& id, const unsigned int& n_species_);
    ~Particle();

    Particle operator+= (const Particle& rhs)& {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xvec[iv] += rhs.xvec[iv];
        }
        return *this;
    }

    friend Particle operator+ (Particle lhs, const Particle& rhs) {
        lhs += rhs;
        return lhs;
    }

    Particle operator-= (const Particle& rhs)& {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xvec[iv] -= rhs.xvec[iv];
        }
        return *this;
    }

    friend Particle operator- (Particle lhs, const Particle& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Particle operator*= (const double& rhs)& {
        for (int iv = 0; iv < n_state_variables; iv++) {
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

    Particle operator/= (const double& rhs)& {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xvec[iv] /= rhs;
        }
        return *this;
    }

    friend Particle operator/ (Particle lhs, const double& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend Particle operator/ (const double& lhs, Particle rhs) {
        rhs /= lhs;
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

    int& getID() {
        return id;
    }

    int& getInjID() {
        return inj_id;
    }

    double& getAge() {
        return age;
    }

    double& getTauRes() {
        return tau_res;
    }

    double& getMass() {
        return mass;
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
        if (setGasState(gas)) {
            return 1.0;
        }
        return gas->density();
    }

    double rho(Cantera::ThermoPhase* gas) {
        if (setGasState(gas)) {
            return 1.0;
        }
        return gas->density();
    }

    double T(std::shared_ptr<Cantera::ThermoPhase> gas) {
        if (setGasState(gas)) {
            return 300.0;
        }
        return gas->temperature();
    }

    double T(Cantera::ThermoPhase* gas) {
        if (setGasState(gas)) {
            return 300.0;
        }
        return gas->temperature();
    }

    double Z(std::shared_ptr<Cantera::ThermoPhase> gas,
             const std::string& comp_fuel,
             const std::string& comp_ox) {
        if (setGasState(gas)) {
            return 0.0;
        }
        return gas->mixtureFraction(comp_fuel, comp_ox);
    }

    double Z(Cantera::ThermoPhase* gas,
             const std::string& comp_fuel,
             const std::string& comp_ox) {
        if (setGasState(gas)) {
            return 0.0;
        }
        return gas->mixtureFraction(comp_fuel, comp_ox);
    }

    void setID(int id_) {
        id = id_;
    }

    void setInjID(int inj_id_) {
        inj_id = inj_id_;
    }

    void setnSpecies(int n_species_) {
        n_species = n_species_;
        n_state_variables = n_species + c_offset_Y;
        xvec.resize(n_state_variables);
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
        for (int k = 0; k < n_species; k++) {
            xvec[c_offset_Y + k] = Y[k];
        }
    }

    void setState(double h, const double* Y) {
        seth(h);
        setY(Y);
    }

    void setState(const double* state) {
        for (int iv = 0; iv < n_state_variables; iv++) {
            xvec[iv] = state[iv];
        }
    }

    int setGasState(std::shared_ptr<Cantera::ThermoPhase> gas) {
        gas->setState_PY(P(), Y());
        try {
            gas->setState_HP(h(), P());
        }
        catch (Cantera::CanteraError) {
            std::cout << "WARNING: Particle " << id << " failed to set state." << std::endl;
            return 1;
        }
        return 0;
    }

    int setGasState(Cantera::ThermoPhase* gas) {
        gas->setState_PY(P(), Y());
        try {
            gas->setState_HP(h(), P());
        }
        catch (Cantera::CanteraError) {
            std::cout << "WARNING: Particle " << id << " failed to set state." << std::endl;
            return 1;
        }
        return 0;
    }

    void print(double threshold = 1.0e-14, std::shared_ptr<Cantera::ThermoPhase> gas=nullptr);

    void react(Cantera::ReactorNet* rnet, double dt);

protected:
    int id;
    int inj_id;
    int n_species;
    int n_state_variables;
    double* m_P = nullptr;
    std::vector<double> xvec;
    double age;
    double mass;
    double tau_res;
    int n_recycles;

private:
};
