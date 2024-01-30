#pragma once
#include "Coupling.h"
#include <map>
#include <utility>
#include "Background.h"

namespace couplings{

/// Photopion production interactions are simulated using tables build with SOPHIA 2.0 generator
/// The code for building the tables is located in sophia subdirectory
/// The tables are adopted for binning density (number of bins per decade) and therefore must be built
/// for each density used
/// The tables should be placed in folder tables/sophia2/<n-bins-per-decade>
class GZK : public Coupling, CouplingParameters
{
	class SecDistribution
	{
	public:
		SecDistribution(const std::string& aPath, bool aBinaryMode = false);
		double N(int deltaE, double aK) const;
		double CalculateRate(const CBackgroundIntegral& aBackgr, const IFunction& aTotCs, int deltaE, double aGamma) const;
		inline int MinDeltaE() const { return fGlobalMinDeltaE; }
		inline int MaxDeltaE() const { return fGlobalMaxDeltaE; }
		inline int NbinsPerDecade() const { return fNbinsPerDecade; }
		inline TParticle Primary() const { return fPrimary; }
		inline TParticle Secondary() const { return fSecondary; }
		inline double Kmin() const { return fKmin; }
		inline double Kmax() const { return fKmax; }
	private:
		double Nbin(int deltaE, int kBin) const;
		CAutoDeletePtrArray<vector<double> > fData;
		vector<int> fMinDeltaE;
		TParticle	fPrimary;
		TParticle	fSecondary;
		double		fKmin;
		double		fKmax;
		double		fStep;
		unsigned int  fNbinsPerDecade;
		unsigned long fNtrials;
		int			fGlobalMinDeltaE;
		int			fGlobalMaxDeltaE;
	};
	class SecKern : public IFunction
	{
	public:
		SecKern(const SecDistribution& aSec, const IFunction& aCs):fCurDeltaE(0),fSec(aSec),fCs(aCs){}
		double f(double _x) const
		{
			return fCs(_x)*fSec.N(fCurDeltaE, _x);
		}
		inline void SetDeltaE(int aDeltaE){fCurDeltaE = aDeltaE;}
	private:
		int fCurDeltaE;
		const SecDistribution& fSec;
		const IFunction& fCs;
	};



	class Channel_N_x : public CouplingChannelT<GZK>
	{
	public:
		Channel_N_x(GZK* aCoupling, TParticle aPrim, TParticle aSec): CouplingChannelT<GZK>(aCoupling, aPrim, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:

	};

	class Channel_A_x : public CouplingChannelT<GZK>
	{
	public:
		Channel_A_x(GZK* aCoupling, TParticle aPrim, TParticle aSec) : CouplingChannelT<GZK>(aCoupling, aPrim, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_A_A : public CouplingChannelT<GZK>
	{
	public:
		Channel_A_A(GZK* aCoupling, TParticle aPrim) : CouplingChannelT<GZK>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	GZK(void);
	~GZK();
	static void GetMultGZKcel(double& multN, double& multP, TParticle aPrim);
	static void GetMultGZKsec(double& multN, double& multP, TParticle aPrim);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	virtual void DebugOutput(const char* aFolder) const;
private:
	void init(double aMaxBackgroundE);
	static const TParticle			secondaries[];
	static const TParticle			primaries[];
	CAutoDeletePtrArray<SecDistribution>	fSecDistributions[2];
	SafePtr<CLogScaleLinearFunc>	fSigmas[2];
	std::map<int, int>				fMapP;
	CAutoDeletePtrArray<CMatrix>	fP[2];
	CVector							fR[2];
	CVector							fCEL[2];
	static bool NO_PI_PRODUCTS;//ignore input of photopions to leptons and photons spectra
	static bool NO_PI_NUCLEI;//do not take into account pion production for nuclei
};									 

}//namespace couplings{
