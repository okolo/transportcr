
#include "Units.h"
#include <math.h>
#include <iostream>

Units::Units(double aEunitInMeV, double aOutEunitInMeV):
plank_Mev_s(6.582122020e-22), 	// Plank constant (MeV*s)
plank_ESU(1.05457266e-27),	// Plank constant (erg*s)
Mpc_cm(3.0856775807e24),		//  (sm)  */
lightspeed(2.99792458e10),	// (sm/s) */
MeVinESU(1.60217733e-6),	// 1MeV in Ergs */
YearInSec(31557600),		//number of seconds in 1 year
MeV_in_K(1.16e10),			//1MeV in K (Kelvin)

outEunitMeV(aOutEunitInMeV),				// Output data energy unit in MeV
Eunit(aEunitInMeV),// units of energy used (in MeV)
Eunit3(aEunitInMeV*aEunitInMeV*aEunitInMeV),// Eunit^3
Tunit(plank_Mev_s/Eunit),// time unit in sec.
Lunit(lightspeed*Tunit),// length unit in sm.
Vunit(Lunit*Lunit*Lunit),// Vunit=Lunit^3
sigmaUnit(Lunit*Lunit),//* sigmaUnit=Lunit^2
SpecUnit(Eunit*Tunit*Vunit),// SpecUnit=Eunit*Tunit*Vunit
Bunit(MeVinESU*MeVinESU*pow(plank_ESU*lightspeed,-1.5)*Eunit*Eunit),// Magnetic field internal unit in Gauss (esu)
barn(1e-24/sigmaUnit),//	1 barn (cross-section unit) in internal units
mbarn(1e-3*barn),
gram(lightspeed*lightspeed/MeVinESU/Eunit),
MeV(1./Eunit),
eV(1e-6*MeV),
keV(1e-3*MeV),
GeV(1e3*MeV),
TeV(1e6*MeV),
PeV(1e9*MeV),
EeV(1e12*MeV),
erg(MeV/MeVinESU),
K(MeV/MeV_in_K),
sec(1./Tunit),
cm(1./Lunit),
cm3(1./Vunit),
Mpc(Mpc_cm*cm),
kpc(Mpc*1e-3),
Hz_photon(2.*M_PI/sec),
J(1e7*erg),
W(J/sec)
{
}

const Units units(1.,1e-6);
