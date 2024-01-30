
#include "RadioBackground.h"
#include <math.h>

//static const double defaultE_rad_max=3e-5;
//static const double defaultE_rad_min=5e-13;

DulkRadioBackground::DulkRadioBackground(bool aIncludeGalactic, bool aIsComovingConst):
	includeGalactic(aIncludeGalactic),
	minFreq(aIncludeGalactic?0.00012:0.001),
	maxFreq(1e4),
	IsComovingConst(aIsComovingConst)
{}

const double DulkRadioBackground::eVtoMHz = 2.41803e8;

double DulkRadioBackground::F(double E, double z)
{
	double gamma1 = E*eVtoMHz;//E in MHertz
	if((gamma1<minFreq)||(gamma1>maxFreq))
		return 0.;
	double tau = 5./pow(gamma1,2.1);
	double ex = Exp(-tau);
	double result = (0.670728/pow(gamma1,0.8)*ex); // EG input
	if(includeGalactic)
	{
		result += (1.56925041509434/pow(gamma1,0.52)*(1.-ex)/tau); // Galactic input
	}
	if(!IsComovingConst)
	{
		double dl = 1.+z;
		result /= (dl*dl*dl);//constant physical concentration (not comoving)
	}
	return result;
}

//returns maximal background red shift
double DulkRadioBackground::MaxZ() const
{
	return 1e100;//infinity (evolution with z not implemented, assuming constant in comoving frame)
}

//returns maximal background energy in eV
double DulkRadioBackground::MaxE(double aZmax) const
{
	return maxFreq/eVtoMHz;
}

//returns minimal background energy in eV
double DulkRadioBackground::MinE(double aZmax) const
{
	return minFreq/eVtoMHz;
}

ProtheroeRadioBackground::ProtheroeRadioBackground(bool aUpper, bool aIsComovingConst):
	y(aUpper ? y_up : y_low),
	IsComovingConst(aIsComovingConst)
{
}

const int ProtheroeRadioBackground::dataSize = 20;

//radio spectrum at z=0 : R.J. Protheroe, P.L. Biermann  astro-ph/9605119 sm^{-3} Fig 5b
//TODO: digitize plot with better accuracy
const double ProtheroeRadioBackground::x[]=
	{
	-6.0,-5.8,-5.6,-5.4,-5.2,-5.0,-4.8,-4.6,-4.4,-4.2,-4.0,-3.8,-3.6,-3.4,-3.0,-2.6,-2.4,-2.2,-1.8,-1.0
};
const double ProtheroeRadioBackground::y_low[]=	// lower curve ( pure luminosity evolution only for radio galaxies )
	{
	-23.3,-23.0,-22.3,-21.9,-21.6,-21.25,-21.1,-20.8,-20.65,-20.6,-20.52,-20.5,-20.5,-20.5,-20.55,-20.58,-20.63,-20.7,-20.95,-21.5
};
const double ProtheroeRadioBackground::y_up[]=	// upper curve ( pure luminosity evolution )
                {
	-22.65,-22.2,-21.85,-21.4,-21.1,-20.7,-20.5,-20.2,-20.1,-19.9,-19.85,-19.8,-19.8,-19.85,-19.95,-20.15,-20.25,-20.3,-20.6,-21.3
};

//Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
//(must be multiplied by (1+z)^3 before substituting farther)
double ProtheroeRadioBackground::F(double E, double z)
{
	const double lightspeed=2.99792458e8; /* (m/s) */
	const double plank=6.582122020e-16; /* Plank constant (eV/s)   */
	const double Pi2=2.*M_PI;
	const double eVinJ=1.60217733e-19;  /* 1eV in J */


	double lg_freq=log10(E/Pi2/plank)-9.0;//log10(frequency/GHz)
	double result;
	if(lg_freq<=x[0])
		result=y[0]+(lg_freq-x[0])*(y[2]-y[0])/(x[2]-x[0]);
	else
		if(lg_freq>=x[dataSize-1])
			result=y[dataSize-1]+(lg_freq-x[dataSize-1])*(y[dataSize-1]-y[dataSize-2])/(x[dataSize-1]-x[dataSize-2]);
		else
		{
			int i;
			for(i=1;x[i]<lg_freq;i++);
			result=y[i-1]+(lg_freq-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1]);
		};
	result=pow(10.0,result); // W*m^{-2}*Hz*{-1}*sr*{-1}
	result/=(1e6*lightspeed*eVinJ*plank*0.5);
	if(!IsComovingConst)
	{
		double dl = 1.+z;
		result /= (dl*dl*dl);//constant physical concentration (not comoving)
	}
	return result;
}

//returns maximal background red shift
double ProtheroeRadioBackground::MaxZ() const
{
	return 1e100;//infinity (evolution with z not implemented, assuming constant in comoving frame)
}

//returns maximal background energy in eV
double ProtheroeRadioBackground::MaxE(double aZmax) const
{
	return 3e-5;
}

//returns minimal background energy in eV
double ProtheroeRadioBackground::MinE(double aZmax) const
{
	return 5e-13;
}
