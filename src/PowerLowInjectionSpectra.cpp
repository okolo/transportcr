// TestInjectionSpectra.cpp: implementation of the CTestInjectionSpectra class.
//
//////////////////////////////////////////////////////////////////////

#include "PowerLowInjectionSpectra.h"
#include "Parameters.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPowerLowInjectionSpectra::CPowerLowInjectionSpectra():
iUseCommonAlpha(true),
iGCRCompositionRates(false),
iCommonAlpha(2.),
iAlphas(EEndParticle),
iNormalizationEnergy(1e13)//MeV
{
	Reader()->readBoolPar("UseCommonAlpha",iUseCommonAlpha);
	Reader()->readBoolPar("GCRCompositionRates",iGCRCompositionRates);
	Reader()->readDoublePar("CommonAlpha",iCommonAlpha);
	Reader()->readDoublePar("NormalizationEnergy",iNormalizationEnergy);
	FOR_ALL_REAL_PARTICLES(particle)
	{
		iAlphas[particle]=2.;//default value
		string settingName = "alpha_";
		settingName += ParticleData::getParticleFileName((TParticle)particle);
		Reader()->readDoublePar(settingName.c_str() ,iAlphas[particle]);
	}
}

CPowerLowInjectionSpectra::~CPowerLowInjectionSpectra()
{

}

double CPowerLowInjectionSpectra::Q(TParticle aParticle, int aBinE, double aE, double aZ)
{
	double alpha = iUseCommonAlpha?iCommonAlpha:iAlphas[aParticle];

	return pow(aE/iNormalizationEnergy,-alpha);
}

void CPowerLowInjectionSpectra::init()
{
	if(iGCRCompositionRates)
	{
		int particle, A;
		double FeRate = GCRrates[56]*pow(56., iCommonAlpha - 1.);
		double norm = FeRate>0?1./FeRate:1.;//normalizing to make iron rate 1 if possible
		for(A=2, particle=(int)EA02;particle<(int)EEndNuclei;particle++, A++)
		{
			s_rates[particle] = norm * GCRrates[A] * pow(A, iCommonAlpha - 1.);
		}
	}
}

const double CPowerLowInjectionSpectra::GCRrates[] = 
// Galactic cosmic rays composition was taken from table 1 col. 4 of paper by Duvernois, M. A. and Thayer, M. R., 1996, ApJ, 465, 982
{
	0, //n
	0, //p
	0, //D
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	4.474,//C12
	0,
	0.342,//N14
	0,
	5.263,//O16
	0,
	0,
	0,
	0.58,//Ne20
	0,
	0,
	0.0324,//Na23
	1.08,//Mg24
	0,
	0,
	0.0778,//Al27
	1.0,//Si28
	0,
	0,
	0.0066,//P31
	0.131,//S32
	0,
	0,
	0.0011,//Cl35
	0,
	0,
	0,
	0.0029,//K39
	0.0873,//Ca40 + Ar40
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0.015,//Cr52
	0,
	0,
	0.011,//Mn55
	0.97//Fe56
};
