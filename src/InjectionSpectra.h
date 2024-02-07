// InjectionSpectra.h: interface for the CInjectionSpectra class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(INJECTIONSPECTRA_H__5F5DD020_C191_11D5_885E_444553540001__INCLUDED_)
#define INJECTIONSPECTRA_H__5F5DD020_C191_11D5_885E_444553540001__INCLUDED_

#include "Sigma.h"
#include "Concentrations.h"
#include "Function.h"

enum TZDependence{
	EZDependenceOff = 0,//(1+z)^3
	EZDependencePowerLaw,//(1+z)^(3+m)
	EZDependencePowerTD,//t^{-4+p}

//injection spectrum types below are proportional to (1+z)^m:

	WaxmanBahcall,//if m=0 corresponds to Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
	SFR03,//if m=0 corresponds to Star Formation Rate from astro-ph/0309141
	Weak,//if m=0 corresponds to m=3 to z=1.8 and constant to z=3 going to zero there
	Baseline,//if m=0 corresponds to baseline evolution astro-ph-1512479
	Strong,//if m=0 corresponds to fast evolution astro-ph-1512479
	Engel9,//if m=0 corresponds to Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
	PowerCut,//(1+z)^m from Zmin up to PowerCutZ and const from PowerCutZ to Zmax
	SFR06,//star formation rate from astro-ph/0607087 (lower figure 3)
	SFR_Yuksel,//star formation rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 11)
	GRB_Yuksel,//GRB rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 12)
	AGN_Aharonian,//AGN rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 13) almost equal to AGN_Hasinger44, probably contains a misprint
	AGN_Hasinger42,//AGN rate from http://arxiv.org/pdf/astro-ph/0506118v1.pdf  (formulas 14,15,16) for Lx=42.5
	AGN_Hasinger43,//AGN rate from http://arxiv.org/pdf/astro-ph/0506118v1.pdf  (formulas 14,15,16) for Lx=43.5
	AGN_Hasinger44,//AGN rate from http://arxiv.org/pdf/astro-ph/0506118v1.pdf  (formulas 14,15,16) for Lx=44.5
	AGN_Hasinger45,//AGN rate from http://arxiv.org/pdf/astro-ph/0506118v1.pdf  (formulas 14,15,16) for Lx=45.5
	Starburst1,
	Starburst2,
	Starforming1,
	Starforming2,
	BLLac,
    EZDependenceCustomTable,//user defined table function of z
	SFR_Madau14,
	EZDependenceDecay,
	EZDependenceCustomCode,

	EZDependenceEOF
};

class CInjectionSpectra  : protected Parameters, public TReferencedObj
{
public:
	/**
	Initializes aSource with this spectrum multipied by aRate
	*/
	CInjectionSpectra(bool aEnableCutOff=true);
	virtual ~CInjectionSpectra();

	//The injection spectrum with all cuts, rate multipliers and density evolution factors applied
	//Parameter aE for nuclei is actually E/A
	//[ModifiedSpectrum] = 1/MeV/sec/cm^3
	double ModifiedSpectrum(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ);

	virtual void init(){};
	virtual double MaxZ() const;
	virtual double MaxE(TParticle aParticle) const;//MeV

	virtual void NormalizeIfNeeded();
	virtual double BaseEmax();//MeV

protected:

	//parameter aE for nuclei is actually E/A
	//injection spectrum without cuts
	virtual double Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ)=0;
	virtual double DensityEvolution(double aZ);

	double SourceCutoffFactor(double aE/*MeV*/, TParticle aParticle);
	static double PowerLawZSourceDependence(double aZ);
	static double PowerCutZSourceDependence(double aZ);
	static double TDZSourceDependence(double aZ);


	static double injSpectraLowerCutoff;
	static double injSpectraLowerAbsCutoff;
	static double injSpectraLowerCutoffWidth;
	static double injSpectraHigherCutoff;
	static double injSpectraPowerCut;
    static double injSpectraPowerCutEnergy;
	static double s_rates[EEndAllParticles];
	static double BlackHoleMass;//black hole mass in units of solar mass
private:
	IFunction* CreateDefaultEvolution();
	enum TMaxEMode{
		ECommonMaxE = 0,// maximal energy for all particles is Emax
		EZMaxE,// maximal energy for nuclei particles is Z*Emax (Hillas criterion) and for the rest of particles Emax
		EAMaxE,// maximal energy for nuclei particles is A*Emax and for the rest of particles Emax
		ECustomMaxE,// custom maximal energy
		EBlackHoleOptimistic,//optimistic model from astro-ph/0808.0367
		EBlackHoleRealistic,//realistic model from astro-ph/0808.0367
		EDiffusive, //Emax~(A/Z)^4 - diffusive acceleration with synchrotron losses (shock wave etc.)
		EInductiveSynchr, //A^2/Z^(3/2) - inductive acceleration with synchrotron-dominated losses (jets of powerful active galaxies)
		EInductiveCurv,  //A/Z^(1/4) - inductive acceleration with curvature-dominated losses

		EMaxEModeEOF //must be the last one
	};
	enum TZminMode{
	EZminModeFixed = 0,
	EZminModeBlackHole
	};

	static bool   s_autoNucleiRates;
	
	static TZDependence		s_zDependence;
	static TZminMode		s_ZminMode;
	static double s_PowerCutZ;
	static bool   s_isInitialized;
	static bool	  s_isSaved;
	static double s_mZ;
	static double s_Zmax;
	static double s_Zmin;
	static double s_p;//TD spectra evolution coefficient
	static double s_Norm;
	static int	  s_MaxEMode;
	double		  m_LastZ;
	double		  m_LastDensity;
	bool          m_enableCutOff;
	static double s_ProtonEmaxMeV;
	SafePtr<IFunction>	m_DefaultEvol;

	class TotalFlux : public IFunction{
	public:
		TotalFlux(CInjectionSpectra* aSpectra):m_spectra(aSpectra){};
		double f(double aE) const;
	private:
		CInjectionSpectra*	m_spectra;
	};

	class AdjustableZDependence : public IFunction
	{
	public:
		AdjustableZDependence(double aDep(double)):m_Dep(new CFunc(aDep)){};
		AdjustableZDependence(IFunction* aDep/*takes ownership*/) : m_Dep(aDep){};
		double f(double aZ) const
		{
			return PowerLawZSourceDependence(aZ)*m_Dep->f(aZ);
		}
	private:
		SafePtr<IFunction> m_Dep;
	};

	/// http://arxiv.org/pdf/astro-ph/0506118v1.pdf (formulas (14),(15),(16) with parameters from table 6)
	class AGN_Evolution_Hasinger : public IFunction
	{
	public:
		AGN_Evolution_Hasinger(double aM, double aZc,double aZd, double aK);
		virtual double f(double aZ) const;
	private:
		double fM;
		double fZc;
		double fZd;
		double fK;
		double fConst;
	};

	class Decay_Evolution : public IFunction
	{
	public:
		Decay_Evolution(double aTau, double aZ_0);
		virtual double f(double aZ) const;
		double fTau;
		double fT0;
	};
};

class CFixedEnergyInjection : public CInjectionSpectra{
public:
	CFixedEnergyInjection();
	virtual double BaseEmax();//MeV
	virtual double MaxE(TParticle aParticle) const;//MeV
	virtual double Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ);
private:
	double iEnergy;
};

#endif // !defined(AFX_INJECTIONSPECTRA_H__5F5DD020_C191_11D5_885E_444553540001__INCLUDED_)
