#include <math.h>
#include "DecaySpectrum.h"
#include "Addfunc.h"
#include "Fragmentation.h"
#include "TableInjSpectra.h"
#include "FragmentationBasedInjSpectra.h"
#include "TimeZ.h"
#include "ParticleList.h"
#include "Ranges.h"

const double MuDecaySpectrum::r=0.573108541682942; /*   r=(Mmuon/Mpion)^2   */
const double MuDecaySpectrum::A0=0.94486;
const double MuDecaySpectrum::A2=-2.7892;
const double MuDecaySpectrum::A3=1.2397;
const double MuDecaySpectrum::B0=-2.4126;
const double MuDecaySpectrum::B_1=-2.8951;
const double MuDecaySpectrum::B2=4.3426;
const double MuDecaySpectrum::B3=-1.9300;
const double MuDecaySpectrum::C0=1.1053;
const double MuDecaySpectrum::C2=-4.46883;
const double MuDecaySpectrum::C3=3.71887;
const double MuDecaySpectrum::D0=13.846;
const double MuDecaySpectrum::D_1=5.37053;
const double MuDecaySpectrum::D1=-28.1116;
const double MuDecaySpectrum::D2=20.0558;
const double MuDecaySpectrum::D3=-5.7902;

MuDecaySpectrum::MuDecaySpectrum()
{
	e0r=-1.*B_1;
	e2r=-2.*B2;
	e3r=-3.*B3;
	en0r=-1.*D_1;
	en2r=-2.*D2;
	en3r=-3.*D3;
	e2l=A2/(1./r/r-1.);
	e3l=A3/(1./r/r/r-1.);
	en2l=C2/(1./r/r-1.);
	en3l=C3/(1./r/r/r-1.);
	e0l=e0r-r*r*(e2r-e2l)-r*r*r*(e3r-e3l);
	en0l=en0r-r*r*(en2r-en2l)-r*r*r*(en3r-en3l);
//test	

	mun_norm=1.;
	aen_norm=1.;

	double eNo=TestMeanNo(EElectron);//3.88
	double aenNo=TestMeanNo(ENeutrinoAE);//3.46
	
	mun_norm=1./eNo;
	aen_norm=1./aenNo;
}

double MuDecaySpectrum::e(double _ratio)
{
	double result;
	double r2=_ratio*_ratio;
/*	if(_ratio<r)
		result=e0r+e2r*r2+e3r*r2*_ratio;
	else*/
		result=e0l+e2l*r2+e3l*r2*_ratio;
	//assert(result>=0.);
	if(result<0.)
	{
		LOG_ERROR4("CMuDecaySpectrum::e(", _ratio, ")=", result);
		return 0.;
	}
	return result*mun_norm;
}
	
double MuDecaySpectrum::aen(double _ratio)
{
	double result;
	double r2=_ratio*_ratio;
/*	if(_ratio<r)
		result=en0r+en2r*r2+en3r*r2*_ratio;
	else*/
		result=en0l+en2l*r2+en3l*r2*_ratio;
	if(result<0.)
	{
		LOG_ERROR4("MuDecaySpectrum::aen(", _ratio, ")=", result);
//		scanf("%*c");
		return 0.;
	}
	return result*aen_norm;
}

double DecaySpectrum::TestMeanNo(TParticle secondaryParticle)
{
	const double dx=1e-5;
	double result=0.;
	for(double x=0.;x<=1.;x+=dx)
		result+=N(x,secondaryParticle);
	result*=dx;
	return result;
}

/**
Energy is conserved if result close to 1
*/
double DecaySpectrum::testEnergyConservation()
{
	double result = 0.;
	FOR_ALL_REAL_PARTICLES_INVOLVED(particle)
	{
		const double dx=1e-5;
		double sum=0.;
		for(double x=0.5*dx;x<=1.;x+=dx)
			sum+=N(x,(TParticle)particle)*x;
		sum*=dx;

		result += sum;
	}
	return result;
}

//moved from Pi_prod.cpp
/*8888888888888888888888888   Pi decay spectra functions  888888888888888888888*/

const double PiDecaySpectrum::r=0.573108541682942; /*   r=(Mmuon/Mpion)^2   */
/*here _z==E_lepton/E_pion,
Ni(E_lepton,E_pion)===Ni(_z)/E_pion,
i=0,1,2,3-the decay channels  */

// Pi0 decay Pi-> 2 gamma :
double PiDecaySpectrum::N0(double _z)
{
 return (_z<1.0)?2.0:0.0;
}
// Pi- ->e-  + ... (the same is valid for muon neutrino prod. spectra)
double PiDecaySpectrum::N1(double _z)
{
 const double A0=0.94486;
 const double A2=-2.7892;
 const double A3=1.2397;
 const double B0=-2.4126;
 const double B_1=-2.8951;
 const double B2=4.3426;
 const double B3=-1.9300;
 if(_z>1.||_z<0.) return 0.;
 double result = ((_z<=r)?(A0+A2*_z*_z+A3*_z*_z*_z):(B0+B_1*log(_z)+B2*_z*_z+B3*_z*_z*_z))/(1.-r);
 double mult = 1.00007;// normalization (1 electron per pion) added 14/05/2002, changed 01/07/2005
 //ASSERT(result>=0);
 return (result>0.)?result*mult:0.;
}
// Pi- ->muon antineutrino + ...
double PiDecaySpectrum::N2(double _z)
{
 return (_z<(1-r))?1.0/(1-r):0;
}
// Pi- ->electron antineutrino + ...
double PiDecaySpectrum::N3(double _z)
{
// it seems that original spectrum is incorrect - it contains step at _z=r and total normalization is wrong
// both problems are solved here temporary by introducing two additional coeffitients (multLeft and multTot)
 const double C0=1.1053;
 const double C2=-4.46883;
 const double C3=3.71887;
 const double D0=13.846;
 const double D_1=5.37053;
 const double D1=-28.1116;
 const double D2=20.0558;
 const double D3=-5.7902;
 if(_z>1.||_z<0.) return 0.;
 double multLeft = 0.71929;//multiplier to avoid step at _z=r
 double multTot = 1.20616;//adjusting overall normalization
 double result = ((_z<=r)?(multLeft*(C0+C2*_z*_z+C3*_z*_z*_z)):D0+D_1*log(_z)+D1*_z+D2*_z*_z+D3*_z*_z*_z)/(1-r);
 return result*multTot;
}

double PiDecaySpectrum::e(double _ratio)
{
	return (m_charge<0)?N1(_ratio):0.;
}

double PiDecaySpectrum::ae(double _ratio)
{
	return (m_charge>0)?N1(_ratio):0.;
}

double PiDecaySpectrum::en(double _ratio)
{
	return (m_charge>0)?N3(_ratio):0.;
}

double PiDecaySpectrum::aen(double _ratio)
{
	return (m_charge<0)?N3(_ratio):0.;
}

double PiDecaySpectrum::mn(double _ratio)
{
	if(!m_charge)
		return 0.;
	return (m_charge<0)?N1(_ratio):N2(_ratio);
}

double PiDecaySpectrum::amn(double _ratio)
{
	if(!m_charge)
		return 0.;
	return (m_charge>0)?N1(_ratio):N2(_ratio);
}

double PiDecaySpectrum::ph(double _ratio)
{
	return (m_charge)?0.:N0(_ratio);
}


