/*
 * BlackHoleInjSpectra.h
 *
 *  Created on: Jan 12, 2012
 *      Author: ok
 */

#ifndef BLACKHOLEINJSPECTRA_H_
#define BLACKHOLEINJSPECTRA_H_

#include "InjectionSpectra.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

class CBlackHoleInjSpectra: public CInjectionSpectra {
public:
	CBlackHoleInjSpectra();
	virtual ~CBlackHoleInjSpectra();
protected:
	//parameter aE for nuclei is actually E/A
	//the method is not thread-safe!
	virtual double Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ);

	void NormalizeIfNeeded();
	void init();
private:
	double EvolutionFactor(double aLog10BHMass);

	//take into account that closest source is located at some z>0
	bool VicinityCut(double aConcentrationMpc3);

	//dN/dLogM [Mpc^-3] (currently arbitrary units)
	//black hole mass is measured in units of solar masses:
	//aLog10BHMass = log10(M/M_solar)
	static double BHMassFunction(double aLog10BHMass);

	//proton Luminosity in erg/sec (currently arbitrary units) for given black hole mass
	//black hole mass is measured in units of solar masses:
	//aLog10BHMass = log10(M/M_solar)
	double Luminosity_p(double aLog10BHMass);

	//fraction of luminosity in gamma to luminosity in p for given black hole mass
	//black hole mass is measured in units of solar masses:
	//aLog10BHMass = log10(M/M_solar)
	double GammaLuminosityFraction(double aLog10BHMass);

	//Monochromatic spectrum is assumed for fixed black hole mass
	//black hole mass is measured in units of solar masses:
	//calculate required black hole mass log10(M/M_solar) for given energy of proton
	//[E] = MeV
	double log10Mp(double aE);

	//Calculate characteristic synchrotron radiation energy for given black hole mass
	//aLog10BHMass = log10(M/M_solar)
	//[E] = MeV
	double gammaEcrit(double aLog10BHMass);

	//Inverse function for gammaEcrit()
	double log10MgammaCrit(double aE);

	//d(Log10(M/M_solar)/dE)
	//Monochromatic spectrum is assumed for fixed black hole mass
	//black hole mass is measured in units of solar masses:
	//calculate required black hole mass log10(M/M_solar) for given energy of photon
	//[E] = MeV
	double dLog10M_dEp(double aE);

private:
	static double SynchroIntencity(double aE, void* aParam);
	static double PhotonKern(double aLog10BHMass, void* aParam);
	double Integration_qags (
			gsl_function aFunction,
			double aXmin,
			double aXmax,
			double epsabs,
			double epsrel,
			size_t limit);

	//χ - angle between magnetic field direction and BH rotation axes ~ 1 deg.
	//currently it is assumed to have the same value for all BHs
	double fChi;

	//Luminosity_p = M^fProtonLuminosityPowerLaw
	double fProtonLuminosityPowerLaw;

	bool fUseVicinityCut;

	static const int SynSize;
	static const double SynX[];
	static const double SynY[];
	double fLowerCoef;
	double fGammaNorm;
    gsl_interp_accel *fAcc;
    gsl_spline *fSpline;
    gsl_integration_workspace* fGslQAGintegrator;
    double fCurE;
    double fCurZ;
    static const double Log10MbhMin;
    static const double Log10MbhMax;
    double fEdingtonBcoef;//magnetic field B=B_edington*fEdingtonBcoef*(M/M_solar)^fEdingtonBalpha
    double fEdingtonBalpha;//magnetic field B=_Bedington*fEdingtonBcoef*(M/M_solar)^fEdingtonBalpha
};

#endif /* BLACKHOLEINJSPECTRA_H_ */
