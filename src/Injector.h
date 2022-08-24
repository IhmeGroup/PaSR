#include <iostream>

#include "cantera/base/ctexceptions.h"

class Injector {
public:
    explicit Injector();
    ~Injector();

    double& T() {
        return m_T;
    }

    double* Y() {
        return m_Y;
    }

    double& Y(const int& k) {
        return m_Y[k];
    }

    void setnsp(const int& nsp_) {
        nsp = nsp_;
    }

    void setT(const double& T_) {
        m_T = T_;
    }

    void setY(const double* Y_) {
        for (int k = 0; k < nsp; k++) {
            m_Y[k] = Y_[k];
        }
    }

    void setState(const double& T, const double* Y) {
        setT(T);
        setY(Y);
    }

protected:
    int nsp;
    double m_T;
    double* m_Y;

private:

};