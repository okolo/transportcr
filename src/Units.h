#ifndef UNITS_H_INCLUDED
#define UNITS_H_INCLUDED

class Units
{
public:
	Units(double aEunitInMeV, double aOutEunitInMeV);
	const double plank_Mev_s; /* Plank constant (MeV*s)   */
	const double plank_ESU;/* Plank constant (erg*s)   */
	const double Mpc_cm;     /*  (cm)  */
	const double lightspeed; /* (cm/s) */
	const double MeVinESU;  /* 1MeV in Ergs */
	const double YearInSec;//number of seconds in 1 year
	const double MeV_in_K;//1MeV in K (Kelvin)
	const double outEunitMeV; /* Output data energy unit in MeV*/
	const double Eunit; /* units of energy used (in MeV) */
	const double Eunit3; /* Eunit^3 */
	const double Tunit;   /* time unit in sec.  */
	const double Lunit;   /* length unit in sm. */
	const double Vunit;   /* Vunit=MpcLunit^3   */
	const double sigmaUnit;  /* sigmaUnit=Lunit^2   */
	const double SpecUnit;  /* SpecUnit=Eunit*Tunit*Vunit */
	const double Bunit;  /* Magnetic field internal unit in Gauss (esu) */

	//////   conventional units expressed in internal units

	const double barn;/*	1 barn (cross-section unit) in internal units */
	const double mbarn;/*	1 milibarn (cross-section unit) in internal units */
	const double gram;// 1 gram in internal units
	const double MeV;
	const double eV;
	const double keV;
	const double GeV;
	const double TeV;
	const double PeV;
	const double EeV;
	const double erg;
	const double K;//Kelvin
	const double sec;
	const double cm;
	const double cm3;//cm^3
	const double Mpc;
	const double kpc;
	const double Hz_photon;// energy of 1 Hz EM wave photon (in internal units)
    const double J; //Joule
    const double W; //Watt
};

extern const Units units;

#endif //#ifndef UNITS_H_INCLUDED //end of file
