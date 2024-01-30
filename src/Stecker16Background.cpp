//
// Created by ok on 01/09/2016.
//

#include "Stecker16Background.h"
#include "Units.h"

Stecker16LowerBackground::Stecker16LowerBackground():
        Stecker16Background(DATA_DIR "stecker_ebl16/comoving_enerdens_lo.csv")
{
}

Stecker16UpperBackground::Stecker16UpperBackground():
        Stecker16Background(DATA_DIR "stecker_ebl16/comoving_enerdens_up.csv")
{
}

/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume */
double Stecker16Background::F(double E, double z){
    return fMatrixFunction.f(z,E)/units.Hz_photon*units.erg;
}

