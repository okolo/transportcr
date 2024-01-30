#include "Kneiske0309IROSpectrum.h"
#include "Units.h"

static const char* files[] = {"0","0.2","0.4","0.6","1","2","3","4",0};
Kneiske0309IROSpectrum::Kneiske0309IROSpectrum():
TableBackground("Kneiske0309", files, true)
{

}

//convert aE to internal x scale
double Kneiske0309IROSpectrum::scaleX(double aE/*eV*/, double aZ)
{
	aE *= (1e-6/units.Eunit);//aE in internal units
	double lambda = 2.*Pi/aE;//wavelength in internal units
	lambda *= units.Lunit*1e4;//wavelength in microns
	return log10(lambda);
}

double Kneiske0309IROSpectrum::unscaleX(double aX, double aZ)
{
	double lambda = pow(10.,aX);//wavelength in microns
	lambda /= units.Lunit*1e4;//wavelength in internal units
	return 2.e6*M_PI/lambda*units.Eunit;
}

//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
double Kneiske0309IROSpectrum::unscaleY(double aY, double aE/*eV*/, double aZ)
{//nW m^-2 sr^-1
	aE *= (1e-6/units.Eunit);//aE in internal units
	double y = pow(10.,aY);
	y *= (4.*Pi*1e-9/(units.MeVinESU*1e-7/units.Eunit)*1e-4*units.Lunit*units.Lunit*units.Tunit);
	double result = y/aE/units.Vunit;
	return result;
}

ElmagKneiskeBestFit::ElmagKneiskeBestFit(bool aIsComoving):
		MatrixBackground(DATA_DIR "kneiskeElmagBestFit.dat", aIsComoving, true)
{
}
