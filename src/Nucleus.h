// Nucleus.h: interface for the CNucleus class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NUCLEUS_H__AD516518_9425_4091_88BF_0844D8105233__INCLUDED_)
#define AFX_NUCLEUS_H__AD516518_9425_4091_88BF_0844D8105233__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ParticleData.h"
#include "Function.h"

class CBackgroundIntegral;

struct TGaussCrossSection{
	double kTh;//threshold energy (MeV)
	double k; //middle energy (MeV)
	double sigma; //total sigma integrated
	double deltaMinus; //energy range (MeV)
	double deltaPlus; //energy range (MeV)
};


struct TBranching{
	int		noOfModes;
	const double* nucleiRates;
	double  pRate;
	double	nRate;
};

struct TIsotope{
	int A; //mass
	int Z; //electric charge

	// photodisintegration with single and double nucleon emission
	TGaussCrossSection pd[2];

	// photodisintegration with multiple nucleon emission
	double pmn;

	// branching rates for multiple nucleon emission
	const TBranching* branching;
};

enum TPhotoDisintegrationChanalType{
	ESingle = 0,//order and numbers are important (see for ex. CPhotoDisintegrationCanal::recalculate)
	EDouble,
	EMultiple
};

class CPhotoDisintegrationCanal{
public:
	CPhotoDisintegrationCanal(const TIsotope& aIsotope, TPhotoDisintegrationChanalType aType);
	
	inline double	R(int iGamma) const{return iR[iGamma];};

	//forces recalculation of R next time R(int) is called
	//used to recalculate R when background changes
	void		recalculate(const CBackgroundIntegral* aBackground);
	void		approximate(CPhotoDisintegrationCanal& aCanal1, CPhotoDisintegrationCanal& aCanal2, double aPart);
	int			getMaxDeltaA();
	double		getMeanDeltaA();//needed for energy loss outputs
	double		getRate(int aDeltaA);
	double		getRateP();
	double		getRateN();
	void		getEnergyLoss(CVector& aOutput);

	inline TParticle	getParticle() const {return (TParticle)(EStartNuclei + iIsotope.A - 2);};
	inline TPhotoDisintegrationChanalType type(){return iType;};
	inline const TIsotope& isotope(){return iIsotope;};
	
public:

	bool isValid();

	class TPDSigma:  public IFunction{
	public:
		virtual double minK()=0;
		virtual double maxK()=0;
	};
		

	class TGaussianSigma :  public TPDSigma{
	public:
		TGaussianSigma(const TGaussCrossSection& aCS);
		virtual double f(double _x) const;//sigma devided by TRK sum
		virtual double minK();
		virtual double maxK();

	private:
		double						iNorm;
		const TGaussCrossSection&	iCS;
		double						iMaxK;
	public:
		static const double			iDefaultMaxK;//30MeV
	};

	class TLinearSigma :  public TPDSigma{
	public:
		TLinearSigma(double  aSigmaIntegral);
		virtual double f(double _x) const{return iNorm;};//sigma devided by TRK sum
		virtual double minK();
		virtual double maxK();

	private:
		double						iNorm;
	public:
		static const double			iMaxK;//150MeV
	};
	
	const TIsotope&					iIsotope;
	CVector							iR;
	TPhotoDisintegrationChanalType	iType;
	double							iSumTRK;
	double							iMeanDeltaA;
};

struct TPhotoDisintegrationMapEntry{
	CPhotoDisintegrationCanal*	iCanal;
	double						iRate;
	TPhotoDisintegrationMapEntry(CPhotoDisintegrationCanal*	aCanal,	double	aRate):iCanal(aCanal),iRate(aRate){};
};

typedef CAutoDeletePtrArray<TPhotoDisintegrationMapEntry>   CPDMSecLine;
typedef CAutoDeletePtrArray<CPhotoDisintegrationCanal>		CPDMLine;//canals for single isotop


class CPhotoDisintegrationMap{
public:
	CPhotoDisintegrationMap();
	
	//reset all entries
	//used to recalculate R when background changes
	void recalculate(const CBackgroundIntegral* aBackground);

	void approximate(CPhotoDisintegrationMap& aMap1, CPhotoDisintegrationMap& aMap2, double aPart);

	// get list of canals contributing to nucleus A
	// A=1 proton
	// A=0 neutron
	const CPDMSecLine* getIncome(int A);

	// get all suppression canals for nucleus A
	const CPDMLine*	getOutcome(int A);
	void printValidity(int aA);
//debug
	void printEnergyLoss(int aA, string& aFileName);

private:

	void initSecMap(CPhotoDisintegrationCanal* aCanal);

	typedef CAutoDeletePtrArray<CPDMLine>						CPDMTable;
	typedef CAutoDeletePtrArray<CPDMSecLine>					CPDMSecTable;
	
	CPDMTable fCanals;
	CPDMSecTable fSecondaryCanals;
};

class Medium;

class CNucleus  /*: public CParticle*/{
public:
	//CNucleus(TParticle aType, bool aNoSpectra = false);

	//virtual ~CNucleus();
	//virtual void saveSpectrum(string dir);
	//virtual void readSpectrum(string dir);
	
//  aA - nucleus mass or A=1 for proton and A=0 for neutron
	//static double incomeCoef(int aA, class Concentrations* _conc, CPhotoDisintegrationMap* aMap, int _i, double _part, double _E);
	//inline double PPPmultiplierP(){return iPPPmultiplierP;};
	static double PPPmultiplierR(int aA);
	static int getZ(TParticle aNucleus);
	static int getA(TParticle aNucleus);
	static void getPhotopionMultipliers(int aA, double& aRn, double& aRp, double& aRN);

	//virtual double getFlux(double aE, int aiE=-1);
/*
	
protected:
	int iZ;//number of protons
	int iA;//number of nucleons
	int iN;//number of neutrons
	double iPPPmultiplierR;
	double iPPPmultiplierP;
	double iPhotopionMultiplierRn;
	double iPhotopionMultiplierRp;
	double iPhotopionMultiplierRN;
	CVector* iEnergies;
	*/
};
/*
template <int A> class CNucleus  : public CNucleusBase  
{
public:
	CNucleus();
	virtual ~CNucleus();
};

//include iplementation
#include "Nucleus.inl"
*/
#endif // !defined(AFX_NUCLEUS_H__AD516518_9425_4091_88BF_0844D8105233__INCLUDED_)
