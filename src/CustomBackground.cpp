//
// Created by ok on 10/13/23.
//

#include "CustomBackground.h"
#include <cmath>

CustomBackground::CustomBackground() :
        fEmin(100),
        fEcut(1e6)
{
}

bool CustomBackground::init()
{
    return true;
}

/* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
  (must be multiplied by (1+z)^3 before substituting farther)*/
double CustomBackground::F(double E, double z)
{
    double e = E/fEcut;
    return exp(-e) / e; // E^{-1} * exp (-E/Ecut)
}

//returns maximal background red shift
double CustomBackground::MaxZ() const
{
    return 1000; // ignore dependence on z
}

//returns maximal background energy in eV
double CustomBackground::MaxE(double aZmax) const
{
    return 10.*fEcut; // we have exponential cutoff at fEcut, so 10*fEcut is enough
}

//returns minimal background energy in eV
double CustomBackground::MinE(double aZmax) const
{
    return fEmin;
}
