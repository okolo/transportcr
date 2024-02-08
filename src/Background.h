#ifndef BACKGROUND_H_INCLUDED
#define BACKGROUND_H_INCLUDED
#include "TableFunc.h"
#include "Vector.h"

extern double E_opt_min;//infrared/optic background minimal energy [eV]

class CTableFunction;

class IBackgroundSpectrum : public TReferencedObj{
protected:
    IBackgroundSpectrum():
    is_initialized(false){}
public:
    inline bool init_once(){
        if (is_initialized)
            return true;
        is_initialized = init();
        return is_initialized;
    }
	/* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z) = 0;

	//returns maximal background red shift
	virtual double MaxZ() const = 0;

	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const = 0;

	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const = 0;
protected:
    virtual bool init(){return true;};
private:
    bool is_initialized;
};

class CuttedBackgroundSpectrum : public IBackgroundSpectrum{
    SmartPtr<IBackgroundSpectrum> fBackground;
    double fMinE; // eV
    double fMaxE; // eV
public:
    CuttedBackgroundSpectrum(IBackgroundSpectrum* aSpec, double aMinE=0., double aMaxE=1e100):
            fBackground(aSpec),
            fMinE(aMinE), // eV
            fMaxE(aMaxE)  // eV
            {}
    virtual bool init(){ return fBackground->init_once(); }
    /* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
      (must be multiplied by (1+z)^3 before substituting farther)*/
    virtual double F(double E, double z){
        return (E < fMinE || E>fMaxE) ? 0. : fBackground->F(E, z);
    }

    //returns maximal background red shift
    virtual double MaxZ() const { return fBackground->MaxZ(); }

    //returns maximal background energy in eV
    virtual double MaxE(double aZmax) const{
        double maxE = fBackground->MaxE(aZmax);
        return maxE > fMaxE ? fMaxE : maxE;
    }

    //returns minimal background energy in eV
    virtual double MinE(double aZmax) const{
        double minE = fBackground->MinE(aZmax);
        return minE < fMinE ? fMinE : minE;
    }
};

class BackgroundUtils
{
public:
    // Calculate integral density in cm^-3
	static double CalcIntegralDensity(IBackgroundSpectrum& aBackground, double aZ, double aRelError=1e-3);
    // Calculate integral energy density in eV * cm^-3
	static double CalcIntegralPowerDensity(IBackgroundSpectrum& aBackground, double aZ, double aRelError=1e-3);
};

class HighRedshiftBackgrExtension : public IBackgroundSpectrum
{
public:
	HighRedshiftBackgrExtension(IBackgroundSpectrum* aBackground, double aDeltaZconst, double aPowerLow, double aDeltaZexp, double aZmax=1000);
	virtual bool init();

	virtual double F(double E, double z);
	//returns maximal background red shift
	virtual double MaxZ() const;
	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;
	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;
private:
	SmartPtr<IBackgroundSpectrum> fBackground;
	double fZconst;
	double fPowerLow;
	double fDeltaZexp;
	double fZmax;
	double fInnerZmax;
};

class CModifiedBlackbodySpectrum : public IBackgroundSpectrum{
public:
	CModifiedBlackbodySpectrum(double aTemperature/*[T]=K*/, double aPower, double aDensity, bool aZdependence, double aZmax, double aWidthDecade=8);
	virtual double MaxZ() const { return iZmax; }
	virtual double F(double E, double z);
	virtual double MaxE(double aZmax) const;
	virtual double MinE(double aZmax) const;
	virtual bool init();
protected:
	double iZmax;
	double iPower;
	double iMultiplier;
	double iT;//temperature
	bool   iZdependence;
	double iRightWidthMult;
	double iLeftWidthMult;
	double iDensity;
};

class IBackgroundObserver
{
public:
	virtual void OnBackgroundChange() = 0;
};

class CBackgroundTable :  public TReferencedObj{
public:
	CBackgroundTable(IBackgroundSpectrum* aSpectrum);
	~CBackgroundTable();
	void update(double aZ/*redshift*/);
	inline double F(int aBin) const{return m_F[aBin];};//returns n(K)deltaK for K =  backgroundE()[aBin]
	inline double Fmid(int aBin) const{return m_Fmid[aBin];};//returns n(K)deltaK for K = midBackgroundE()[aBin]
	inline const int nK() const{return m_nK;};
	inline const CBinning& backgroundE() const{return m_backgroundEnergies;};
	inline const CBinning& midBackgroundE() const{return m_midBackgroundEnergies;};
	inline const double Kmax() const{return m_Kmax;};
	inline const double Kmin() const{return m_Kmin;};

	const CLinearFunc* getBackground();
	void saveBackground(const char* dir) const;
	void AddObserver(IBackgroundObserver* aObserver);
private:
	void NotifyObservers();

	CVarArray<IBackgroundObserver*> m_Observers;
	SmartPtr<IBackgroundSpectrum> fSpectrum;
	CVector			m_F;
	CVector			m_Fmid;
	CLinearFunc*	m_Ffit;
	CVector			m_F_all;
	CVector			m_E_all;
	CBinning 		m_backgroundEnergies;
	CBinning 		m_midBackgroundEnergies;
	int				m_nK;
	double 			m_Kmax;
	double 			m_Kmin;
	double			m_curZ;
};

class CompoundBackground : public IBackgroundSpectrum{
public:
    CompoundBackground(double aEmin=0, double aEmax=1e100);
    CompoundBackground(CompoundBackground* aCopyFrom);
	void addComponent(IBackgroundSpectrum* aComponent, double aWeight=1.);
    void replacePart(IBackgroundSpectrum* aComponent, double aWeight=1.);

	virtual bool init();
	virtual double F(double E, double z);
	virtual double MaxZ() const;
	virtual double MaxE(double aZmax) const;
	virtual double MinE(double aZmax) const;

protected:
	CSmartPointerArray<IBackgroundSpectrum>	m_components;
    SmartPtr<IBackgroundSpectrum>           m_overriden_backgr;
    double                                  m_overriden_backgr_weight;
	std::vector<double>						m_weights;
	double                                  m_Emin;
	double                                  m_Emax;
};

ostream& operator<<(ostream& target, const CBackgroundTable& toWrite);

class CBackgroundIntegral : protected IBackgroundObserver{ //public CBackgroundObserversCollection
public:
	CBackgroundIntegral(CBackgroundTable& aBackground);
	virtual ~CBackgroundIntegral();

/************************************************************************/
/* Calculate suppression coefficient based on the cross section in the
	rest frame of the CR particle
	aGamma - gamma factor of the particle in lab frame
	*/
	double calculateR(double aGamma, IFunction* aSigma, double minK, double maxK, int aAc=1000) const;

	/*
	  Calculate suppression coefficient based on the cross section in the
	  rest frame of the CR particle
	  aGamma - gamma factor of the particle in lab frame
	  using GSL integration API
	*/
	double calculateR(double aGamma, IFunction* aSigma, double minK, double maxK, double aRelError, double epsabs=0, int key=GSL_INTEG_GAUSS15) const;

	/*
	  Calculate suppression coefficient based on the cross section in the
	  rest frame of the CR particle
	  aGamma - gamma factor of the particle in lab frame
	  using GSL integration API on log scale
	*/
	double calculateRlogsc(double aGamma, IFunction* aSigma, double minK, double maxK, double aRelError, double epsabs=0, int key=GSL_INTEG_GAUSS15) const;

	/*
	Integrate[n(b)/b^2,{b,aE,Infinity}]
	*/
	double integral(double aE) const;

	///Return total rate kern multiplied by (2.*aGamma^2)
	inline IFunction* CreateTotalRateKern(IFunction* aSigma, double aGamma) const
	{
		return new TBackgroundFunctional(*aSigma,*fTable,aGamma);
	}

	/*
	Integrate[Integrate[n(b)/b^2,{b,x,Infinity}],{x,aE,Infinity}]
	Used for DPP secondaries calculation (in constant sigma approximation)
	*/
	//double integralDpp(double aE) const;
	inline double operator ()(double aE){return integral(aE);};
	inline double Kmin() const {return (*fK)[0];}
	inline double Kmax() const {return (*fK)[fK->length()-1];}
protected:
	virtual void OnBackgroundChange();
	void reset();
	class TBackgroundFunctional : public IFunction{
	public:
		TBackgroundFunctional(const IFunction& aSigma, CTableFunction& aBackgroundIntegral, double aGamma):
		  fBackgroundIntegral(aBackgroundIntegral),fGamma(aGamma),fSigma(aSigma){};
		virtual double f(double _x) const;
	private:
		CTableFunction&		fBackgroundIntegral;
		double				fGamma;
		const				IFunction& fSigma;
	};
	class TBackgroundFunctionalLogsc : public IFunction{
		public:
		TBackgroundFunctionalLogsc(const IFunction& aSigma, CTableFunction& aBackgroundIntegral, double aGamma):
			  fBackgroundIntegral(aBackgroundIntegral),fGamma(aGamma),fSigma(aSigma){};
			virtual double f(double _x) const;
		private:
			CTableFunction&		fBackgroundIntegral;
			double				fGamma;
			const				IFunction& fSigma;
	};

	//the tabulated integral
	CTableFunction*			fTable;
	

	double					fTotIntegral;
	CVector*				fK;
	CVector*				fI;
	CBackgroundTable&		fBackground;
};


#endif //end of header
