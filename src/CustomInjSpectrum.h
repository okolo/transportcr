//
// Created by ok on 5/7/23.
//

#ifndef PROPAGATION_CUSTOMINJSPECTRUM_H
#define PROPAGATION_CUSTOMINJSPECTRUM_H

#include "InjectionSpectra.h"

class CCustomInjSpectrum : public CInjectionSpectra
{
public:
    CCustomInjSpectrum();
    //parameter aE for nuclei is actually E/A
    //injection spectrum without cuts
    virtual double Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ);
    virtual void init();
};

inline double CustomEvolution(double aZ){
    // https://www.overleaf.com/project/6453d74335a88241cb2f84b3
    const double delta = 1.92;
    const double a = 0.0157;
    const double b = 0.118;
    const double u = 3.23;
    const double k = 4.66;
    return pow(1.+aZ,delta)*(a+b*aZ)/(1.+pow(aZ/u, k))/a; //norm was choosen to have 1 at z=0
}

#endif //PROPAGATION_CUSTOMINJSPECTRUM_H
