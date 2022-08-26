#include <iostream>

#include "cantera/base/ctexceptions.h"

class Injector {
public:
    explicit Injector();
    explicit Injector(const int& index, const unsigned int& nsp_);
    ~Injector();

    double& h() {
        return m_h;
    }

    double& T() {
        throw Cantera::NotImplementedError("Injector::setT");
    }

    std::vector<double>& Y() {
        return m_Y;
    }

    double& Y(const int& k) {
        return m_Y[k];
    }

    void setIndex(const int& index_) {
        index = index_;
    }

    void setnsp(const unsigned int& nsp_) {
        nsp = nsp_;
        m_Y.resize(nsp);
    }

    double getFlow() {
        return flow;
    }

    void setFlow(const double& flow_) {
        flow = flow_;
    }

    void seth(const double& h_) {
        m_h = h_;
    }

    void setT(const double& T_) {
        throw Cantera::NotImplementedError("Injector::setT");
    }

    void setY(const double* Y_) {
        for (int k = 0; k < nsp; k++) {
            m_Y[k] = Y_[k];
        }
    }

    void setState(const double& h, const double* Y) {
        seth(h);
        setY(Y);
    }

    void print(double threshold = 1.0e-14);

protected:
    int index;
    unsigned int nsp;
    double flow;
    double m_h;
    std::vector<double> m_Y;

private:

};