/*
 * BlackHoleInjSpectra.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: ok
 */

#include "BlackHoleInjSpectra.h"
#include "Parameters.h"
#include "TimeZ.h"
#include "Units.h"

CBlackHoleInjSpectra::CBlackHoleInjSpectra():
fChi(M_PI/180.),//1 deg.
fProtonLuminosityPowerLaw(7./8.),// two models considered by Sergey Troitsky et al. 3/8 and 7/8
fUseVicinityCut(true),
fLowerCoef(0),
fGammaNorm(0),
fAcc(0),
fSpline(0),
fGslQAGintegrator(0),
fCurE(0),
fCurZ(0),
fEdingtonBcoef(1.),
fEdingtonBalpha(0.)
{
	double chiDegrees = fChi*180./M_PI;
	Reader()->readDoublePar("NeronovBH_chi",chiDegrees);
	fChi=chiDegrees*M_PI/180.;
	Reader()->readDoublePar("NeronovBH_LuminosityPowerLaw",fProtonLuminosityPowerLaw);
	Reader()->readBoolPar("NeronovBH_UseVicinityCut",fUseVicinityCut);
	Reader()->readDoublePar("EdingtonBcoef",fEdingtonBcoef);
	Reader()->readDoublePar("EdingtonBalpha",fEdingtonBalpha);
}

CBlackHoleInjSpectra::~CBlackHoleInjSpectra() {
	if(fGslQAGintegrator)
		gsl_integration_workspace_free (fGslQAGintegrator);
	if(fSpline)
		gsl_spline_free (fSpline);
	if(fAcc)
		gsl_interp_accel_free (fAcc);
}

void CBlackHoleInjSpectra::NormalizeIfNeeded()
{
	bool useVicinityCut = fUseVicinityCut;
	fUseVicinityCut = false;
	CInjectionSpectra::NormalizeIfNeeded();
	fUseVicinityCut = useVicinityCut;
}

double CBlackHoleInjSpectra::Q(TParticle aParticle, int aBinE, double aE, double aZ)
{
	fCurZ = aZ;
	if(aParticle==EProton)
	{
		double log10Mbh = log10Mp(aE);
		double mf = BHMassFunction(log10Mbh);
		double evolFactor = EvolutionFactor(log10Mbh);
		if(evolFactor>0 && VicinityCut(mf))
			return mf*Luminosity_p(log10Mbh)/aE*dLog10M_dEp(aE)*evolFactor;
	}
	else if(aParticle==EPhoton)
	{
		double EcritMin = aE/SynX[SynSize-1];
		double log10MbhMin = log10MgammaCrit(EcritMin);
		if(log10MbhMin>=Log10MbhMax)
			return 0.;
		if(log10MbhMin<Log10MbhMin)
			log10MbhMin=Log10MbhMin;
		fCurE = aE;
	    gsl_function kern;
	    kern.function = PhotonKern;
	    kern.params = this;
		double result = Integration_qags (kern, log10MbhMin, Log10MbhMax, 0, 1e-4, 1000);
		return result;
	}
	return 0.;
}

double CBlackHoleInjSpectra::EvolutionFactor(double aLog10BHMass)
{
	const double a = 4.30079781470395;
	const double b = 7.18476422945172;
	const double c = 198203.649118613;
	const double d = 398504.359573236;
	const double e = -148078;
	if(fCurZ>3.)
		return 0.;
	double result = (c+d*fCurZ+e*fCurZ*fCurZ)/pow(fCurZ+a,b);
	ASSERT_VALID_NO(result);
	return result;
}

bool CBlackHoleInjSpectra::VicinityCut(double aConcentrationMpc3)
{
	if(aConcentrationMpc3==0.)
		return false;
	if(!fUseVicinityCut)
		return true;
	double d = CTimeZ::z2d(fCurZ)*units.Lunit/units.Mpc_cm;
	double l=pow(aConcentrationMpc3,-0.333333333333333);
	return d>l;
}

const double CBlackHoleInjSpectra::Log10MbhMin = 5.;
const double CBlackHoleInjSpectra::Log10MbhMax = 11.;

//dN/dLogM [Mpc^-3] (currently arbitrary units)
//black hole mass is measured in units of solar masses
//aLogBHMass = log10(M/M_solar)
double CBlackHoleInjSpectra::BHMassFunction(double aLog10BHMass)
{// taken from fig. 5 of Mon. Not. R. Astron. Soc. 354, 1020–1030 (2004)
 // parameterized using gnuplot fit command
 // f(x)=a*x**4+b*x**3+c*x**2+d*x+e
 // fit f(x) './bh.dat' using 1:2 via a, b, c, d, e
 // range shown on the plot: 6 < aLogBHMass < 10
	if(aLog10BHMass<Log10MbhMin||aLog10BHMass>Log10MbhMax)
		return 0.;//cut unphysical range
	const double coef[] = {-79.6331,44.4838,-9.55544,0.915862,-0.0332728};
	double sum = coef[0];
	double powerX=aLog10BHMass;
	for(int power=1; power<5; power++, powerX*=aLog10BHMass)
		sum += (coef[power]*powerX);
	return pow(10., sum);
}

//proton Luminosity in erg/sec (currently arbitrary units) for given particle type and black hole mass
//black hole mass is measured in units of solar masses:
//aLog10BHMass = log10(M/M_solar)
double CBlackHoleInjSpectra::Luminosity_p(double aLog10BHMass)
{
	return pow(10., aLog10BHMass*fProtonLuminosityPowerLaw);
}

//fraction of luminosity in gamma to luminosity in p for given black hole mass
//black hole mass is measured in units of solar masses:
//aLog10BHMass = log10(M/M_solar)
double CBlackHoleInjSpectra::GammaLuminosityFraction(double aLog10BHMass)
{
	return 31.7*sqrt(fChi*sqrt(fEdingtonBcoef*fEdingtonBcoef*fEdingtonBcoef))*pow(10.,(0.125+0.75*fEdingtonBalpha)*aLog10BHMass-1.125);
}

//Monochromatic spectrum is assumed for fixed black hole mass
//black hole mass is measured in units of solar masses:
//calculate required black hole mass log10(M/M_solar) for given energy of proton
//[E] = MeV
double CBlackHoleInjSpectra::log10Mp(double aE)
{
	return (3.375+log10(aE*sqrt(fChi/sqrt(fEdingtonBcoef))/2.79e13))/(0.375+0.25*fEdingtonBalpha);
}

//Monochromatic spectrum is assumed for fixed black hole mass
//black hole mass is measured in units of solar masses:
//calculate required black hole mass log10(M/M_solar) for given energy of photon
//[E] = MeV
double CBlackHoleInjSpectra::log10MgammaCrit(double aE)
{
	return (1.125+log10(aE/4.4e6/pow(fEdingtonBcoef,0.75)))/(0.125+0.75*fEdingtonBalpha);
}

double CBlackHoleInjSpectra::gammaEcrit(double aLog10BHMass)
{
    //double fEdingtonBcoef;//magnetic field B=B_edington*fEdingtonBcoef*(M/M_solar)^fEdingtonBalpha
    //double fEdingtonBalpha;//magnetic field B=_Bedington*fEdingtonBcoef*(M/M_solar)^fEdingtonBalpha
	return 4.4e6*pow(fEdingtonBcoef,0.75)*pow(10, -1.125 + (0.125+0.75*fEdingtonBalpha)*aLog10BHMass);
}

//d(Log10(M/M_solar)/dE)
//Monochromatic spectrum is assumed for fixed black hole mass
//black hole mass is measured in units of solar masses:
//calculate required black hole mass log10(M/M_solar) for given energy of photon
//[E] = MeV
double CBlackHoleInjSpectra::dLog10M_dEp(double aE)
{
	return 0.434294481903252/(0.375+0.25*fEdingtonBalpha)/aE;// 1/log(10)/(3/8+fEdingtonBalpha/4)/E
}

void CBlackHoleInjSpectra::init()
{
    fAcc = gsl_interp_accel_alloc ();
    fSpline = gsl_spline_alloc (gsl_interp_linear, SynSize);
    gsl_spline_init (fSpline, SynX, SynY, SynSize);
    fLowerCoef = SynY[0]/pow(SynX[0],0.3333333333333333);
    fGammaNorm = 1.;
    gsl_function phi;
    phi.function = SynchroIntencity;
    phi.params = this;
    double norm = Integration_qags (phi, 0., SynX[SynSize-1], 0, 1e-4, 1000);
    fGammaNorm = 1./norm;//normalize to unit luminosity
}

double CBlackHoleInjSpectra::PhotonKern(double aLog10BHMass, void* aParam)
{
	CBlackHoleInjSpectra* self = (CBlackHoleInjSpectra*)aParam;
	double Ecrit = self->gammaEcrit(aLog10BHMass);
	double x = self->fCurE/Ecrit;
	double mf = self->BHMassFunction(aLog10BHMass);
	double evolFactor = self->EvolutionFactor(aLog10BHMass);

	if(evolFactor>0 && self->VicinityCut(mf))
	{
		double phi = SynchroIntencity(x, self);
		return mf*self->Luminosity_p(aLog10BHMass)*self->GammaLuminosityFraction(aLog10BHMass)*evolFactor*phi/Ecrit/self->fCurE;
	}
	return 0.;
}

double CBlackHoleInjSpectra::SynchroIntencity(double aX, void* aParam)
{
	CBlackHoleInjSpectra* self = (CBlackHoleInjSpectra*)aParam;
	if(aX<SynX[0])
		return self->fGammaNorm*self->fGammaNorm*self->fLowerCoef*pow(aX,0.3333333333333333);
	if(aX>SynX[SynSize-1])
		return 0.;

	return self->fGammaNorm*gsl_spline_eval (self->fSpline, aX, self->fAcc);
}

double CBlackHoleInjSpectra::Integration_qags (
		gsl_function aFunction,
		double aXmin,
		double aXmax,
		double epsabs,
		double epsrel,
		size_t limit)
{
	if(fGslQAGintegrator==0)
		fGslQAGintegrator = gsl_integration_workspace_alloc (limit);
	else if(fGslQAGintegrator->limit < limit)
	{
		gsl_integration_workspace_free (fGslQAGintegrator);
		fGslQAGintegrator = gsl_integration_workspace_alloc (limit);
	}
	double result, abserr;
	int failed = gsl_integration_qags (&aFunction, aXmin, aXmax, epsabs, epsrel, limit, fGslQAGintegrator, &result, &abserr);

	if(failed)
	{
		ASSERT(0);
		ThrowError("Integration failed with code " + ToString(failed));
	}
	return result;
}

const int CBlackHoleInjSpectra::SynSize = 71;

const double CBlackHoleInjSpectra::SynX[] = {
		1.e-6,
		1.2589254117941661e-6,
		1.584893192461114e-6,
		1.9952623149688787e-6,
		2.5118864315095823e-6,
		3.162277660168379e-6,
		3.981071705534969e-6,
		5.011872336272725e-6,
		6.30957344480193e-6,
		7.943282347242822e-6,
		0.00001,
		0.000012589254117941661,
		0.00001584893192461114,
		0.000019952623149688786,
		0.000025118864315095822,
		0.000031622776601683795,
		0.000039810717055349695,
		0.00005011872336272725,
		0.00006309573444801929,
		0.00007943282347242822,
		0.0001,
		0.00012589254117941674,
		0.00015848931924611142,
		0.0001995262314968881,
		0.0002511886431509582,
		0.00031622776601683794,
		0.00039810717055349735,
		0.0005011872336272725,
		0.0006309573444801936,
		0.0007943282347242822,
		0.001,
		0.0012589254117941675,
		0.001584893192461114,
		0.0019952623149688807,
		0.002511886431509582,
		0.0031622776601683794,
		0.003981071705534973,
		0.005011872336272725,
		0.006309573444801936,
		0.00794328234724282,
		0.01,
		0.012589254117941687,
		0.01584893192461114,
		0.01995262314968879,
		0.025118864315095822,
		0.03162277660168379,
		0.039810717055349776,
		0.05011872336272725,
		0.06309573444801943,
		0.07943282347242822,
		0.1,
		0.12589254117941687,
		0.15848931924611143,
		0.19952623149688828,
		0.25118864315095824,
		0.31622776601683794,
		0.39810717055349776,
		0.5011872336272725,
		0.6309573444801942,
		0.7943282347242822,
		1.,
		1.2589254117941688,
		1.584893192461114,
		1.9952623149688828,
		2.5118864315095824,
		3.1622776601683795,
		3.9810717055349776,
		5.011872336272725,
		6.309573444801943,
		7.943282347242821,
		10.
};

const double CBlackHoleInjSpectra::SynY[] = {
		0.021493468615984588,
		0.023207788621204248,
		0.02505878464722503,
		0.02705733826087169,
		0.02921519346705981,
		0.03154502404942224,
		0.034060505904881806,
		0.03677639467264665,
		0.039708608956723526,
		0.04287431943298675,
		0.04629204411487714,
		0.049981750023085476,
		0.05396496146046541,
		0.058264875029305795,
		0.06290648143810648,
		0.06791669402170472,
		0.07332448373257718,
		0.0791610201404643,
		0.0854598176870618,
		0.09225688606334322,
		0.09959088308506678,
		0.10750326780682562,
		0.11603845079832883,
		0.1252439374602,
		0.13517045891963583,
		0.14587208334223525,
		0.1574062983288749,
		0.16983405231461526,
		0.1832197393982824,
		0.1976311076160318,
		0.21313906509145036,
		0.22981735145104695,
		0.24774203301733316,
		0.2669907691242932,
		0.2876417828792938,
		0.3097724521368027,
		0.33345741453137806,
		0.3587660531757033,
		0.3857591959680557,
		0.4144848201945227,
		0.44497250411421074,
		0.47722630762570956,
		0.5112156948339709,
		0.546864033893457,
		0.584034128561529,
		0.6225101618682812,
		0.6619753857128425,
		0.7019849087434067,
		0.7419330846263715,
		0.7810153949129041,
		0.8181855348728537,
		0.8521099241362534,
		0.8811244783460038,
		0.9032027203897995,
		0.9159507527614681,
		0.9166536000581921,
		0.9024084135965875,
		0.8703902498786638,
		0.8182984877410949,
		0.7450125112233047,
		0.6514228153553645,
		0.5412769858116798,
		0.42169259427242983,
		0.3028163346598772,
		0.19614237756062256,
		0.11150081123878199,
		0.053736564846518155,
		0.021017612757718383,
		0.006314173497486309,
		0.0013594602058937969,
		0.00019223826521752585
};
