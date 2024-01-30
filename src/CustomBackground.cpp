//
// Created by ok on 10/13/23.
//

#include "CustomBackground.h"
#include "Units.h"
#include <cmath>

CustomBackground::CustomBackground() :
        fEmin(100),
        fEcut(1e6),
        fNormN(1.)
{

}

bool CustomBackground::init()
{
    double norm_luminosity = 1e44 * units.erg / units.sec;  // visible luminosity of the source in internal units
    double radius = 1e13 * units.cm;  // radius of the region filled with the radiation
    double nE = fNormN / units.cm3 * fEcut * log(fEcut/fEmin); // energy density in internal units assuming dn/dE ~ 1/E^2
    double luminosity = 2*M_PI*radius*radius*nE; // luminosity in internal units before normalization
    fNormN *= norm_luminosity/luminosity;  // norm factor to achieve visible luminosity
    fNormN /= 3.7753932507E-10; // temporary fix
    fNormN  *= 0.00044;
    return true;
}

/* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
  (must be multiplied by (1+z)^3 before substituting farther)*/
double CustomBackground::F(double E, double z)
{
    double e = E/fEcut;
    return fNormN * exp(-e) / e; // E^{-1} * exp (-E/Ecut)
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
