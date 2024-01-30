// TestInjectionSpectra.h: interface for the CTestInjectionSpectra class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TESTINJECTIONSPECTRA_H__BFE086A2_0968_4F1B_B024_B88984DCB12F__INCLUDED_)
#define AFX_TESTINJECTIONSPECTRA_H__BFE086A2_0968_4F1B_B024_B88984DCB12F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "InjectionSpectra.h"

class CPowerLowInjectionSpectra  :  public CInjectionSpectra
{
public:
	CPowerLowInjectionSpectra();
	virtual ~CPowerLowInjectionSpectra();

	virtual double Q(TParticle aParticle, int aBinE, double aE, double aZ);
	virtual void init();
private:
	static const double GCRrates[];
	bool	iUseCommonAlpha;
	bool	iGCRCompositionRates;
	double	iCommonAlpha;
	CVector iAlphas;
	double	iNormalizationEnergy;//MeV
};

#endif // !defined(AFX_TESTINJECTIONSPECTRA_H__BFE086A2_0968_4F1B_B024_B88984DCB12F__INCLUDED_)
