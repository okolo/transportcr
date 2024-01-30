//
// Created by ok on 29/05/18.
//

#include "Franceschini17EBL.h"

static const char* files[] = {"0", "0.2", "0.4", "0.6", "0.8", "1", "1.2", "1.4", "1.8", "2", 0};
Franceschini17EBL::Franceschini17EBL():
        TableBackground("Franceschini17", files, true)
{

}

//convert aE to internal x scale
double Franceschini17EBL::scaleX(double aE/*eV*/, double aZ)
{
    return log10(aE);
}

//convert internal x scale to energy in eV
double Franceschini17EBL::unscaleX(double aX, double aZ)
{
    return pow(10.,aX);
}

//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
double Franceschini17EBL::unscaleY(double aY, double aE/*eV*/, double aZ)
{
    double dV = (1.+aZ);
    dV *= (dV*dV);
    return pow(10., aY)/dV; // tables are given for physical concentration in cm^-3
}

