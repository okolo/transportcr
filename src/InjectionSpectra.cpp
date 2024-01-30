// InjectionSpectra.cpp: implementation of the CInjectionSpectra class.
//
//////////////////////////////////////////////////////////////////////

#include "InjectionSpectra.h"
#include "CustomInjSpectrum.h"
#include "Parameters.h"
#include "Ranges.h"
#include "ParticleList.h"
#include <math.h>
#include "TimeZ.h"
#include "Nucleus.h"
#include "Function.h"
#include "Units.h"
#include "PropagEngine.h"
#include <fstream>
#include "TableFunc.h"
#include "TableReader.h"
#include "BLLac.h"

inline double weakEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
//m=3 to z=1.8 and constant to z=3 going to zero there
	return (aZ<1.8)?pow(1.+aZ,3.):((aZ<3.)?21.952:0.);
}

//baseline evolution astro-ph-0512479
inline double baselineEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
//m=4 to z=1 and then constant to z=6 going to zero there.
	return (aZ<1.3)?pow(1.+aZ,3.1):((aZ<6)?13.2238005912547:0.);
}

//fast evolution astro-ph-0512479
inline double strongEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
//m=4 to z=1 and then constant to z=6 going to zero there.
	return (aZ<1.0)?pow(1.+aZ,4.):((aZ<6)?16.:0.);
}

//AGN
//Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
inline double WaxmanBahcallEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	return (aZ<1.9)?pow(1.+aZ,3.):(24.389*((aZ<2.7)?1.:Exp((2.7-aZ)/2.7)));
}

//Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
inline double Engel9Evolution(double aZ){//corresponds to comoving frame only if m=0!!!
	return ((aZ<1.9)?pow(1.+aZ,4.):70.7281);
}

//Star Formation Rate from astro-ph/0309141
inline double SfrEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	return (aZ<1.2)?pow(1.+aZ, 3.5):40.68052976*pow(1.+aZ, -1.2);//40.68052976==(1+1.2)^(3.5+1.2)
}

//star formation rate from astro-ph/0607087 (lower figure 3)
inline double SfrNewEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	if(aZ<1.) 
		return pow(1.+aZ, 4.65);
	else if(aZ<3.)
		return 7.7274906314*pow(1.+aZ, 1.7);
	else if(aZ<6.)
		return 7912.950406552335*pow(1.+aZ, -3.3);
	else if(aZ<8.)
		return 606880463687.06*pow(1.+aZ, -12.63);
	else
		return 0.;
}

//star formation rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 11)
inline double SFR_YukselEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	if(aZ<1.)
		return pow(1.+aZ, 3.4);
	else if(aZ<4.)
		return 10.5560632861832*pow((1.+aZ)/2., -0.3);
	else
		return 8.01899573803639*pow((1.+aZ)/5., -3.5);
}

//GRB rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 12)
inline double GRB_YukselEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	if(aZ<1.)
		return pow(1.+aZ, 4.8);
	else if(aZ<4.)
		return 27.857618025476*pow((1.+aZ)/2., 1.1);
	else
		return 76.3269641062938*pow((1.+aZ)/5., -2.1);
}

//Type II supernova rate from astro-ph/0601463
inline double SNII(double aZ){
	const double a=0.0118;
	const double b=0.08;
	const double c=3.3;
	const double d=5.2;
//const double h=0.7;
	return (a+b*aZ)/(1.+pow(aZ/c,d));
}

//example 1 from T. A. Thompson, E. Quataert, E. Waxman and A. Loeb, astro-ph/0608699
inline double frac_starburst1(double aZ)
{
	return aZ<1.?(0.9*aZ+0.1) : 1.;
}

//example 2 from T. A. Thompson, E. Quataert, E. Waxman and A. Loeb, astro-ph/0608699
inline double frac_starburst2(double aZ)
{
	if(aZ<1.)
	{
		double z1 = aZ+1.;
		return 0.1*z1*z1*z1;
	}
	else
		return 0.8;
}

inline double StarburstGalaxies1(double aZ){
	return frac_starburst1(aZ)*SNII(aZ);
}

inline double StarburstGalaxies2(double aZ){
	return frac_starburst2(aZ)*SNII(aZ);
}

inline double StarFormingGalaxies1(double aZ){
	return (1.-frac_starburst1(aZ))*SNII(aZ);
}

inline double StarFormingGalaxies2(double aZ){
	return (1.-frac_starburst2(aZ))*SNII(aZ);
}

inline double SFR_Madau14rate(double aZ){
	double z1=1.+aZ;
	return 1.00257377875529*pow(z1,2.7)/(1.+pow(z1/2.9,5.6));//norm was choosen to have 1 at z=0
}

/// http://arxiv.org/pdf/astro-ph/0506118v1.pdf (formulas (14),(15),(16) with parameters from table 6)
CInjectionSpectra::AGN_Evolution_Hasinger::AGN_Evolution_Hasinger(double aM, double aZc,double aZd, double aK):
fM(aM),fZc(aZc),fZd(aZd),fK(aK)
{
	fConst = pow(1.+aZc, aM);
}

double CInjectionSpectra::AGN_Evolution_Hasinger::f(double aZ) const
{
	double mult = PowerLawZSourceDependence(aZ);
	if(aZ<fZc)
		return mult*pow(1.+aZ, fM);
	else if(aZ<fZd)
		return mult*fConst;
	else
		return mult*fConst*pow(10.,fK*(aZ-fZd));
}

CInjectionSpectra::Decay_Evolution::Decay_Evolution(double aTau, double aZ_0):
        fTau(aTau),
        fT0(CTimeZ::z2t(aZ_0))
{
}

double CInjectionSpectra::Decay_Evolution::f(double aZ) const
{
    double dt = CTimeZ::z2t(aZ)-fT0;
    return Exp(-dt/fTau);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CInjectionSpectra::CInjectionSpectra(bool aEnableCutOff):
m_LastZ(-1.),
m_LastDensity(-1.),
m_enableCutOff(aEnableCutOff)
{
	if (s_isInitialized) {
		return;
	}

	FOR_ALL_REAL_PARTICLES(particle){
		const char* name = ParticleData::getParticleFileName((TParticle)particle);
		string parName = name;
		parName += "Rate";
		s_rates[particle] = 0.;//default value is 0
		Reader()->readDoublePar(parName.c_str(),s_rates[particle]);
	}
	Reader()->readBoolPar("AutoNucleiRates", s_autoNucleiRates);
	if(s_autoNucleiRates)
	{
		s_rates[EProton] = 0;
		s_rates[ENeutron] = 0;
		TParticle heaviestNucleus = EEndAllParticles;
		FOR_ALL_NUCLEI_INVOLVED(particle)
		{
			s_rates[particle] = 0;
			heaviestNucleus = particle;
		}
		if(heaviestNucleus != EEndAllParticles)
			s_rates[heaviestNucleus] = 1;
		else if(CParticleList::Instance()->IsEnabled(EProton))
		{
			s_rates[EProton] = 1;
		}
	}

	READ_DOUBLE_SETTING(injSpectraLowerAbsCutoff);
	READ_DOUBLE_SETTING(injSpectraLowerCutoff);
	READ_DOUBLE_SETTING(injSpectraLowerCutoffWidth);
	READ_DOUBLE_SETTING(injSpectraHigherCutoff);
    READ_DOUBLE_SETTING(injSpectraPowerCut);
    READ_DOUBLE_SETTING(injSpectraPowerCutEnergy); // eV


	Reader()->readSwitchPar("EmaxMode", s_MaxEMode);
    Reader()->readDoublePar("BlackHoleMass",BlackHoleMass);
	s_MaxEMode = s_MaxEMode%EMaxEModeEOF;

	switch(s_MaxEMode)
	{
		case EBlackHoleOptimistic:
		{
			double x = BlackHoleMass*pow(10.,-8.45);
			s_ProtonEmaxMeV = 5.46e13/units.Eunit*pow(x, 0.375);
			break;
		}
		case EBlackHoleRealistic:
		{
			double x = BlackHoleMass*pow(10.,-8.45);
			s_ProtonEmaxMeV = 1.85e13/units.Eunit*pow(x, 0.2975);
			break;
		}
		default:
		{
			double protonEmaxEv = 1e21;
			Reader()->readDoublePar("Emax",protonEmaxEv);
			s_ProtonEmaxMeV = protonEmaxEv*1e-6;
		}
	}

	int s_zDependenceInt = s_zDependence;
	Reader()->readSwitchPar("z-dependence", (int&)s_zDependenceInt);
	if(s_zDependenceInt<0 || s_zDependenceInt >= EZDependenceEOF)
		throw "invalid z-dependence parameter";
	s_zDependence = (TZDependence)s_zDependenceInt;

	Reader()->readDoublePar("m_z",s_mZ);
	Reader()->readDoublePar("Zmax",s_Zmax);
	Reader()->readDoublePar("Zmin",s_Zmin);
	Reader()->readDoublePar("PowerCutZ",s_PowerCutZ);
	Reader()->readDoublePar("p",s_p);

	int ZminMode = s_ZminMode;
	Reader()->readSwitchPar("ZminMode", ZminMode);
	s_ZminMode = (TZminMode)ZminMode;

	if(s_ZminMode == EZminModeBlackHole)
	{
		double x = BlackHoleMass*pow(10.,-8.45);
		s_Zmin = 0.0012*pow(x,-0.23)*exp((x-1.)/3.);//the formula is written for H=70km/s/Mpc, todo: introduce explicit dependence on Hubble constant
	}

	s_isInitialized = true;
}

CInjectionSpectra::~CInjectionSpectra()
{

}

double CInjectionSpectra::s_rates[EEndAllParticles];
bool   CInjectionSpectra::s_autoNucleiRates = false;
TZDependence CInjectionSpectra::s_zDependence = EZDependencePowerLaw;
CInjectionSpectra::TZminMode CInjectionSpectra::s_ZminMode = EZminModeFixed;
bool   CInjectionSpectra::s_isInitialized = false;
bool   CInjectionSpectra::s_isSaved = false;
double CInjectionSpectra::s_mZ = 0.;
double CInjectionSpectra::s_Zmax = 3.;
double CInjectionSpectra::s_Zmin = 0.;
double CInjectionSpectra::s_p = 1;
double CInjectionSpectra::injSpectraLowerCutoff = 0.;
double CInjectionSpectra::injSpectraLowerAbsCutoff = 0.;
double CInjectionSpectra::injSpectraLowerCutoffWidth = 3.5;
double CInjectionSpectra::injSpectraHigherCutoff = 2.;
double CInjectionSpectra::injSpectraPowerCut = 0.3;
double CInjectionSpectra::injSpectraPowerCutEnergy = 1e23;
int CInjectionSpectra::s_MaxEMode = CInjectionSpectra::ECommonMaxE;
double CInjectionSpectra::s_PowerCutZ = 2;
double CInjectionSpectra::BlackHoleMass = 1e6;


double CInjectionSpectra::ModifiedSpectrum(TParticle aParticle, int aBinE, double aEnergyInMeV, double z)
{
	double rate = s_rates[aParticle]*s_Norm;
	ASSERT_VALID_NO(rate);
	if (rate == 0.) return 0.;

	double E = aEnergyInMeV*ParticleData::GetEnergyScaleFactor(aParticle);
	if(E<ParticleData::getParticleMass(aParticle)*units.Eunit)
		return 0.;
	if (m_enableCutOff) {
		rate *= SourceCutoffFactor(E, aParticle);
		ASSERT_VALID_NO(rate);
		if (rate == 0.) return 0.;
	}

	if(z!=m_LastZ)//not thread-safe!
	{
		m_LastDensity = DensityEvolution(z);
		m_LastZ = z;
	}
	ASSERT_VALID_NO(m_LastDensity);
	rate *= m_LastDensity;
	if (rate == 0.) return 0.;

	double result = rate*Q(aParticle, aBinE, E, z);
	ASSERT_VALID_NO(result);
	return result;
}

double CInjectionSpectra::DensityEvolution(double z)
{
	if(m_DefaultEvol.isNull())
	{
		m_DefaultEvol = CreateDefaultEvolution();
#ifdef _DEBUG
        {
            std::cerr << "\n\n\n--------------Evolution-----------\n\n";
            for(double oz=0; oz<Ranges().Zmax(); oz+=0.1){
                std::cerr << oz << "\t" << m_DefaultEvol->f(oz) << "\n";
            }
            std::cerr << "\n\n-------------------------------------\n\n";
        };
#endif
	}
	return m_DefaultEvol->f(z);
}

IFunction* CInjectionSpectra::CreateDefaultEvolution()
{
	switch (s_zDependence) 
	{
	case EZDependencePowerLaw:
		return new CFunc(PowerLawZSourceDependence);
	case EZDependencePowerTD:
		return new CFunc(TDZSourceDependence);
	case EZDependenceOff:
		return new CFunc(unitFunc);
	case WaxmanBahcall://if m=0 corresponds to Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
		return new AdjustableZDependence(WaxmanBahcallEvolution);
	case SFR03://if m=0 corresponds to Star Formation Rate from astro-ph/0309141
		return new AdjustableZDependence(SfrEvolution);
	case SFR06://if m=0 corresponds to star formation rate from astro-ph/0607087 (lower figure 3)
		return new AdjustableZDependence(SfrNewEvolution);
	case Weak://if corresponds to (1+z)^(m+3) to z=1.8 and (1+z)^m to z=3 going to zero there
		return new AdjustableZDependence(weakEvolution);
	case Baseline://if m=0 corresponds to baseline evolution astro-ph-1512479
		return new AdjustableZDependence(baselineEvolution);
	case Strong://if m=0 corresponds to fast evolution astro-ph-1512479
		return new AdjustableZDependence(strongEvolution);
	case Engel9://if m=0 corresponds to Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
		return new AdjustableZDependence(Engel9Evolution);
	case SFR_Yuksel:
		return new AdjustableZDependence(SFR_YukselEvolution);
	case GRB_Yuksel:
		return new AdjustableZDependence(GRB_YukselEvolution);
	case AGN_Aharonian:
		///AGN rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 13)
		///equal to AGN_HasingerEvol44_5 at aZ<2.7 but with fK=-1 (probably a misprint in 1103.3574v2)
		return new AGN_Evolution_Hasinger(5.,1.7,2.7,-1.);
	case AGN_Hasinger42://Log(Lx)=42.5
		return new AGN_Evolution_Hasinger(4.,0.7,0.7,-0.32);
	case AGN_Hasinger43://Log(Lx)=43.5
		return new AGN_Evolution_Hasinger(3.4,1.2,1.2,-0.32);
	case AGN_Hasinger44://Log(Lx)=44.5 also used in Ahlers at. al arxiv/0902.3993 formula (10)
		return new AGN_Evolution_Hasinger(5.,1.7,2.7,-0.43);
	case AGN_Hasinger45://Log(Lx)=45.5
		return new AGN_Evolution_Hasinger (7.1,1.7,2.7,-0.43);
	case PowerCut://(1+z)^m from Zmin up to PowerCutZ and const from PowerCutZ to Zmax
		return new CFunc(PowerCutZSourceDependence);
	//TODO:  write comments in corresponding entries of xsw switch section when paper is ready
	case Starburst1:
		return new AdjustableZDependence(StarburstGalaxies1);
	case Starburst2:
		return new AdjustableZDependence(StarburstGalaxies2);
	case Starforming1:
		return new AdjustableZDependence(StarFormingGalaxies1);
	case Starforming2:
		return new AdjustableZDependence(StarFormingGalaxies2);
	case EZDependenceCustomCode:
        return new AdjustableZDependence(CustomEvolution);
	case BLLac:
		{
			//SafePtr<CLogScaleLinearFunc> evol = new CLogScaleLinearFunc(new CTableReader(DATA_DIR "evolution/bllac",2));
			//evol->SetExtension(ExtLinear);
			return new AdjustableZDependence(new BLLacLDDE());
		}
    case EZDependenceCustomTable:
        {
            std::string customEvol = "custom_evol";
            Reader()->readStringPar("CustomEvolFile", customEvol);
            if(!fileExists(customEvol))
                ThrowError("Invalid CustomEvolFile parameter: file not found");
            CopyFileOrDir(customEvol, plt_local_dir);
            //using log scale interpolation in (1+z)
            return new AdjustableZDependence(new ArgShiftFunc(new CLogScaleLinearFunc((new CTableReader(customEvol, 2))->ShiftCol(0, 1.)), 1.));
        }
    case SFR_Madau14:
		{
			return new AdjustableZDependence(SFR_Madau14rate);
		}
    case EZDependenceDecay:
    {
        double SourceDecayT = 1e15;
        Reader()->readDoublePar("SourceDecayT", SourceDecayT);
        return new Decay_Evolution(SourceDecayT * units.sec, s_Zmax);
    }

	default:
		break;
	}
	throw "invalid value of z-dependence switch";
}

double CInjectionSpectra::s_Norm = -1.;
double CInjectionSpectra::s_ProtonEmaxMeV = 0.;

class LogEnergyDensityKern : public IFunction
{
public:
	LogEnergyDensityKern(IFunction& aSpectrum):
	fSpectrum(aSpectrum){}
	double f(double aLogE) const
	{
		double E = exp(aLogE);
		return fSpectrum(E) * E * E;
	}
private:
	IFunction& fSpectrum;
};

class LogDensityKern : public IFunction
{
public:
	LogDensityKern(IFunction& aSpectrum):
			fSpectrum(aSpectrum){}
	double f(double aLogE) const
	{
		double E = exp(aLogE);
		return fSpectrum(E) * E;
	}
private:
	IFunction& fSpectrum;
};

double CInjectionSpectra::TotalFlux::f(double aE_MeV) const
{
	double sum = 0.;
	FOR_ALL_PARTICLES_INVOLVED(particle)
	{
		double scalingFactor = ParticleData::GetEnergyScaleFactor(particle);
		sum+=m_spectra->ModifiedSpectrum(particle, -1, aE_MeV/scalingFactor, 0.);
	}
	return sum;
}

void CInjectionSpectra::NormalizeIfNeeded()
{
	if(s_Norm > 0)
		return;//source was already normalized, skip normalization
	s_Norm = 1.;
	if(s_MaxEMode == EBlackHoleOptimistic || s_MaxEMode == EBlackHoleRealistic || s_ZminMode == EZminModeBlackHole)
	{
		LOG_WARNING("skipping spectrum normalization");
		return;
	}
	NormMode normMode = NormalizeEnergyDensity;
	Reader()->readSwitchPar("SourceNormMode", &normMode, EndNormMode);
	if(normMode == NormalizeAbsolute) {
		Reader()->readDoublePar("SourceAbsNorm", s_Norm);
		if(s_Norm<0)
			ThrowError("Invalid parameter for absolute normalization: SourceAbsNorm=" + ToString(s_Norm));
	}
	double EmaxMeV = Ranges().midE()[Ranges().nE()-1]/units.MeV;
	double SourceIntegrDensity = 2.9379989438e73;// 1/cm3/[yr]
	double SourceIntegrEnergyDensity = 4.5e44;// erg/Mpc3/[yr] (astro-ph/0101216v2)
	double DensityIntegrEmin = 0;
	double DensityIntegrEmax = 0;

	READ_DOUBLE_SETTING(SourceIntegrDensity);
	READ_DOUBLE_SETTING(SourceIntegrEnergyDensity);
	READ_DOUBLE_SETTING(DensityIntegrEmin);
	READ_DOUBLE_SETTING(DensityIntegrEmax);

	double E1MeV = DensityIntegrEmin >0 ? DensityIntegrEmin *1e-6 : Ranges().midE()[0]/units.MeV;
	double E2MeV = DensityIntegrEmax >0 ? DensityIntegrEmax *1e-6 : EmaxMeV;

	if(E1MeV > EmaxMeV)
	{
		LOG_WARNING("skipping spectrum norm calculation (Emax is too low)");
		return;	
	}

	TotalFlux spec(this);
	double zMin = CInjectionSpectra::s_Zmin;
	CInjectionSpectra::s_Zmin = 0;

	SafePtr<IFunction> kern = new LogDensityKern(spec);
	SafePtr<IFunction> kernE = new LogEnergyDensityKern(spec);
	double density = FuncUtils::integrate(log(E1MeV), log(E2MeV), kern, &gUnitFunc, 10000);
	ASSERT_VALID_NO(density);
	double densityE = FuncUtils::integrate(log(E1MeV), log(E2MeV), kernE, &gUnitFunc, 10000);
	ASSERT_VALID_NO(densityE);

	double multDensity = 1./units.Mpc_cm/units.Mpc_cm/units.Mpc_cm;
	if(CPropagEngine::GetSourceType() == PointSource)
	{
		//for moment source N_ini=Q*dt, where dt = SourceTable::DefaultMomentSourceDt in internal units
		multDensity /= (SourceTable::DefaultMomentSourceDt*units.Tunit); //multiplier for [SourceIntegr(Energy)Density] = (erg)/Mpc^3
	}
	else
	{
		multDensity /= units.YearInSec;//multiplier for [SourceIntegr(Energy)Density] = (erg)*Mpc^-3/year
	}
	double multDensityE = multDensity*units.erg/units.MeV;
	if(normMode == NormalizeAbsolute)
	{
		SourceIntegrDensity = density/multDensity;
		SourceIntegrEnergyDensity = densityE/multDensityE;
	}
	else {
		if(density==0.){
			WARN("spectrum normalization failed (zero flux). Using norm 1");
		}
		else {
			if (normMode == NormalizeEnergyDensity) {
				s_Norm = multDensityE * SourceIntegrEnergyDensity / densityE;
				SourceIntegrDensity = s_Norm * density / multDensity;
			}
			else if (normMode == NormalizeDensity) {
				s_Norm = multDensity * SourceIntegrDensity / density;
				SourceIntegrEnergyDensity = s_Norm * densityE / multDensityE;
			}
			else
				ThrowError("Invalid NormMode parameter: " + ToString(normMode));
		}
	}
	ASSERT_VALID_NO(s_Norm);
	ASSERT_VALID_NO(SourceIntegrDensity);
	ASSERT_VALID_NO(SourceIntegrEnergyDensity);
	IParWriter* writer = IParWriter::Instance();
	if(writer)
	{
		if(normMode != NormalizeAbsolute)
			writer->writePar("SourceAbsNorm", s_Norm);
		if(normMode != NormalizeDensity)
			writer->writePar("SourceIntegrDensity", SourceIntegrDensity);
		if(normMode != NormalizeEnergyDensity)
			writer->writePar("SourceIntegrEnergyDensity", SourceIntegrEnergyDensity);
	}
	CInjectionSpectra::s_Zmin = zMin;
}

double CInjectionSpectra::BaseEmax()
{
	switch(s_MaxEMode){
		case EDiffusive://factor pow(A/Z,4)
			return s_ProtonEmaxMeV*8.;// for A=2 Emax=16EmaxP
		case EInductiveSynchr: //A^2/Z^(3/2)
			return s_ProtonEmaxMeV*2.;// for A=2 Emax=4EmaxP
		default:
			return s_ProtonEmaxMeV;
	}
}

double CInjectionSpectra::SourceCutoffFactor(double aE/*MeV*/, TParticle aParticle){
	double maxE = MaxE(aParticle);

	if (aE>maxE) return 0.;

	double result = 1.;

	if(injSpectraHigherCutoff<1)
		result *= Exp(-aE/(injSpectraHigherCutoff*maxE));

	double power_cut = aE*1e6 / injSpectraPowerCutEnergy;
    if(power_cut > 1.)
        result *= pow(power_cut, -injSpectraPowerCut);
	
	if(injSpectraLowerAbsCutoff>0.)
		result *= Exp(-1e-6*injSpectraLowerAbsCutoff/(aE*units.Eunit));
	else
	{
		double R = injSpectraLowerCutoff*maxE/aE;
		if (R>1) {//below cutoff energy
			if(injSpectraLowerCutoffWidth <=0.)
				return 0.;//sharp cutoff
			result *= pow(R,-log10(R)/injSpectraLowerCutoffWidth/injSpectraLowerCutoffWidth);
		}
	}
	return result;
}

double CInjectionSpectra::MaxE(TParticle aParticle) const
{
	double maxE = s_ProtonEmaxMeV;

	if (aParticle>=EStartNuclei&&aParticle<EEndNuclei){
		double A = CNucleus::getA(aParticle);
		double Z = CNucleus::getZ(aParticle);
		switch(s_MaxEMode){
			//case ECommonMaxE:// maximal energy for all particles is Emax
			case EZMaxE:// maximal energy for nuclei particles is Z*Emax and for the rest of particles Emax
				maxE *= Z;
				break;
			case EAMaxE:// maximal energy for nuclei particles is A*Emax and for the rest of particles Emax
				maxE *= A;
				break;
			case ECustomMaxE:// custom maximal energy
				//maxE = MaxE(aParticle);
				ThrowError("CustomMaxE mode is not supported by selected InjectionSpectrum implementation");
				break;
			case EBlackHoleOptimistic:
			case EBlackHoleRealistic:
				maxE *= (A * pow(Z, -0.25));
				break;
			case EDiffusive:
			{
				double ratio = A/Z;
				maxE *= (ratio*ratio*ratio*ratio);
				break;
			}
			case EInductiveSynchr: //A^2/Z^(3/2)
				maxE *= (A*A/Z/sqrt(Z));
				break;
			case EInductiveCurv:  //A/Z^(1/4)
				maxE *= (A/pow(Z,0.25));
				break;
		}
	}
	return maxE;
}

double CInjectionSpectra::MaxZ() const
{
	return s_Zmax;
}

double CInjectionSpectra::PowerLawZSourceDependence(double aZ)
{
	if((aZ<s_Zmin)||(aZ>s_Zmax))
		return 0.0;
	return pow(1.0+aZ, s_mZ);
}

double CInjectionSpectra::PowerCutZSourceDependence(double aZ)
{
	if((aZ<s_Zmin)||(aZ>s_Zmax))
		return 0.0;
	if(aZ>s_PowerCutZ)
		aZ = s_PowerCutZ;
	return pow(1.0+aZ, s_mZ);
}

double CInjectionSpectra::TDZSourceDependence(double aZ)
{
	double dl = 1.+aZ;
	double CosmT = CTimeZ::z2t(aZ)/CTimeZ::getAgeOfUniverse(); /*cosmic time - CosmT==t(sec)/ageOfUniverse   */
	const double A=2.0e-44;     /* dNx/dCosmT=A*CosmT^(-4+p)   */  //for p=1
	return (A*pow(CosmT,-4.0+s_p)/(dl*dl*dl));//converting to comoving conc.
}

double CFixedEnergyInjection::BaseEmax()
{
	return iEnergy;
}

double CFixedEnergyInjection::MaxE(TParticle aParticle) const
{
	return iEnergy;
}

CFixedEnergyInjection::CFixedEnergyInjection():
CInjectionSpectra(false),
iEnergy(0)//MeV
{
	Reader()->readDoublePar("DeltaSourceEnergy",iEnergy);//If 0 or less Emax is used
	if(iEnergy<=0)
		iEnergy = CInjectionSpectra::BaseEmax();
}

double CFixedEnergyInjection::Q(TParticle aParticle, int aBinE, double aE, double aZ)
{
	if(aE<iEnergy*BC().ss_2)
		return 0.;
	if(aE>iEnergy*BC().ss2)
		return 0.;
	return 1.;
}
