// Ranges.h: interface for the CRanges class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RANGES_H__5035C553_CE88_11D5_B789_E1E4F5AD3C9A__INCLUDED_)
#define AFX_RANGES_H__5035C553_CE88_11D5_B789_E1E4F5AD3C9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "main.h"
#include "Vector.h"
#include "TParticle.h"

struct TBinningConstants{
	TBinningConstants(int aNoOfBinsPerDecade);
	int NoOfBinsPerDecade;/* number of bins per decade */
	double s; //same as NoOfBinsPerDecade but stored in double

/* constants which depend on number of bins per decade */
	double ss1; // 10^(1/s)
	double ss2; // sqrt(ss1)
	double ss_1;// 1/ss1
	double ss_2;// 1/ss2
	double ss3; // ss1^3
	double sF;  // ss2-ss_2
	double sLn; // log(ss1)
};

class CBinning : public CVector{
public:
	void Create(double aInitialVal, int aNoOfBins, const TBinningConstants& aBc);
	inline int getBin(double aX) const {
		int result = findX(aX);
		if (result >= 0) return result;
		if (aX<=(*this)[0]) return aX>m_minVal?0:-1;
		if (aX>=(*this)[m_length-1]) return aX<m_maxVal?m_length-1:-1;
		return -1;
	};
private:
	double m_minVal;
	double m_maxVal;
};

class CRanges  
{
public:
	CRanges(int aAccuracy, double aApproxEmin, double aApproxEmax, 
		double aZmax, int aAccuracyZ, bool aCELredshift, double aMinInteractionEnergyMeV,
		double aMinNucleiInteractionEnergyMeV);
	virtual ~CRanges();
	
	inline const CBinning& E() const{return m_energies;};
	inline const CBinning& midE() const{return m_midEnergies;};
	inline const CBinning& nucleonGamma() const{return m_nucleonGammas;};
	inline const CBinning& midNucleonGamma() const{return m_midNucleonGammas;};
	//inline const CBinning& backgroundE() const{return m_backgroundEnergies;};
	//inline const CBinning& midBackgroundE() const{return m_midBackgroundEnergies;};
	inline const CVector& Z() const {return m_Z;}

	inline const double Emax() const{return m_Emax;};
	inline const double Emin() const{return m_Emin;};
	//inline const double Kmax() const{return m_Kmax;};
	//inline const double Kmin() const{return m_Kmin;};
	inline const double Zmax() const{return m_Zmax;};
	inline const double emin(int aTimeStep) const{return m_energies[aTimeStep];};
	inline const TBinningConstants& getBinningConstants() const{return m_bc;};

	inline const int nE() const{return m_nE;};
	//inline const int nK() const{return m_nK;};
	inline const int nMinInteraction(TParticle p=EStartLightParticle) const
	{
	    return p < EStartNuclei ? m_MinInteractionBin : m_MinInteractionBinN;
	};
    inline const int nMinInteractionN() const{return m_MinInteractionBinN;};
	inline static const CRanges& Static(){return *s_defaultRanges;};
	// returns previous default ranges
	static CRanges* SetDefault(CRanges* aRanges);

	inline int Accuracy() const { return m_accuracy;};
	inline int AccuracyZ() const { return m_accuracyZ;};
	
	/// Facilitate creating binning with exactly iAccuracy * s_minNumberOfBinsPerDecade points per decade
	/// in the way to ensure that for any integer n E=10^n eV fits exactly into the middle of a bin
	/// returns number of bins
	int adjustRanges(double& aMin, double& aMax) const;

protected:

	CBinning m_energies;
	CBinning m_midEnergies;
	CBinning m_nucleonGammas;
	CBinning m_midNucleonGammas;


	CVector m_Z;

	double m_Emin;
	double m_Emax;

	TBinningConstants m_bc;

	double m_Zmax;

	int m_nE;//number of CR energy bins

	static const int s_minNumberOfBinsPerDecade;
	static CRanges*  s_defaultRanges;

	int m_accuracy;
	int m_accuracyZ;
	int m_MinInteractionBin;
	int m_MinInteractionBinN;

friend class CTableFunction;
};

//extern const double meanNucleonMass; /* mean nuclon mass (MeV) */

/**
Shotcut functions and definitions
*/

inline const CRanges& Ranges(){
	return CRanges::Static();
}

inline const TBinningConstants& BC(){
	return CRanges::Static().getBinningConstants();
}


#define DECLARE_E_SHORTCUT \
	const double Emin = Ranges().Emin();\
	const double Emax = Ranges().Emax();\
	const int nn = Ranges().nE();

#define DECLARE_K_SHORTCUT \
	const double bmin = Ranges().Kmin();\
	const double bmax = Ranges().Kmax();\
	const int mm = Ranges().nK();


#define DECLARE_BINNING_CONST_SHORTCUT \
	const double ss1 = BC().ss1;\
	const double ss2 = BC().ss2;\
	const double ss_1 = BC().ss_1;\
	const double ss_2 = BC().ss_2;\
	const double ss3 = BC().ss3;\
	const double sF = BC().sF;\
	const double s = BC().s;
	


#endif // !defined(AFX_RANGES_H__5035C553_CE88_11D5_B789_E1E4F5AD3C9A__INCLUDED_)
