#if !defined(PRODSPEC_H_INCLUDED_)
#define PRODSPEC_H_INCLUDED_

#include "ParticleData.h"
#include <assert.h>

class DecaySpectrum  
{
public:

	virtual double e(double){return 0.0;};
	virtual double ae(double){return 0.0;};
	virtual double en(double){return 0.0;};
	virtual double aen(double){return 0.0;};
	virtual double mn(double){return 0.0;};
	virtual double amn(double){return 0.0;};
	virtual double tn(double){return 0.0;};
	virtual double atn(double){return 0.0;};
	virtual double ph(double){return 0.0;};
	virtual double p(double){return 0.0;};
	virtual double n(double){return 0.0;};
	virtual double N(double _ratio,TParticle secondaryParticle)
	{
		switch(secondaryParticle)
		{
			case EElectron:
				return e(_ratio);
			case EPositron:
				return ae(_ratio);
			case EPhoton:
				return ph(_ratio);
			case ENeutrinoE:
				return en(_ratio);
			case ENeutrinoM:
				return mn(_ratio);
			case ENeutrinoT:
				return tn(_ratio);
			case ENeutrinoAE:
				return aen(_ratio);
			case ENeutrinoAM:
				return amn(_ratio);
			case ENeutrinoAT:
				return atn(_ratio);
			case ENeutron:
				return n(_ratio);
			case EProton:
				return p(_ratio);
		default:
		  return 0.;
		};
		return 0.0;
	};
	double N(double _E_primary,double _E_secondary,TParticle secondaryParticle)
	{
		double ratio=_E_secondary/_E_primary;
		if(ratio>1.)
			return 0.;
		return N(ratio,secondaryParticle)/_E_primary;
	};
	
/**
Energy conservation test:
energy is conserved if result close to 1
*/
	double testEnergyConservation();
	double TestMeanNo(TParticle secondaryParticle);
};

class MuDecaySpectrum : public DecaySpectrum
{
//this spectrum is dirived from pion decay spectrum used in this program,
//the last one is presented in the paper of S. Lee astro-ph/9604098
protected:
static const double r; /*   r=(Mmuon/Mpion)^2   */
static const double A0;
static const double A2;
static const double A3;
static const double B0;
static const double B_1;
static const double B2;
static const double B3;
static const double C0;
static const double C2;
static const double C3;
static const double D0;
static const double D_1;
static const double D1;
static const double D2;
static const double D3;
double e0l;
double e2l;
double e3l;
double en0l;
double en2l;
double en3l;
double e0r;
double e2r;
double e3r;
double en0r;
double en2r;
double en3r;

double mun_norm;
double aen_norm;

public:	
	MuDecaySpectrum();
	virtual double e(double _ratio);
	
	virtual double aen(double _ratio);
	
	virtual double mn(double _ratio){return MuDecaySpectrum::e(_ratio);};
};

class TauDecaySpectrum : public DecaySpectrum
{
public:
	virtual double tn(double _ratio){return 1.0;};//temporary
};

class PiDecaySpectrum : public DecaySpectrum
{
public:
	PiDecaySpectrum(int charge/* +1, -1 or 0 */):m_charge(charge){};

	virtual double e(double _ratio);
	virtual double ae(double _ratio);
	virtual double en(double _ratio);
	virtual double aen(double _ratio);
	virtual double mn(double _ratio);
	virtual double amn(double _ratio);
	virtual double ph(double _ratio);
protected:
	static double N0(double _z);
	static double N1(double _z);
	static double N2(double _z);
	static double N3(double _z);
	static const double r;
	int m_charge;/* +1, -1 or 0 */
};

#endif
