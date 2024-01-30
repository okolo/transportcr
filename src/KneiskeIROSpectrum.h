// KneiskeIROSpectrum.h: interface for the CKneiskeIROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_KNEISKEIROSPECTRUM_H__66C50CA8_C861_41BC_A180_A2B3916CC96F__INCLUDED_)
#define AFX_KNEISKEIROSPECTRUM_H__66C50CA8_C861_41BC_A180_A2B3916CC96F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TableBackground.h"
#include "TableFunc.h"

// IR/O spectrum from Tanja Kneiske (kneiske@astro.uni-wuerzburg.de) et al. 
// see http://www.astro.uni-wuerzburg.de/theorie/

class CKneiskeIROSpectrum : public CTableWithHeaderBackgr  
{
public:
	CKneiskeIROSpectrum();
	virtual ~CKneiskeIROSpectrum();
// from CIROSpectrum
	double F(double E, double z);

	double MinE(double aZmax) const;
	double MaxE(double aZmax) const;
	double MaxZ() const;
private:
// from CDataReader
	ECompleteStatus readHeaderLine(const char* theString);
	bool readDataLine(const char* theString);
	bool testData();
	void processData();

// data
	int m_accuracyE;
	int m_accuracyZ;
	int m_curE;
	CVector m_zArray;// redshifts
	CVector m_eArray;// energies
	CMatrix m_nArray;// photon concentrations
	CAutoDeletePtrArray<CTableFunction> m_fE; /// F(E,z_i)
};

#endif // !defined(AFX_KNEISKEIROSPECTRUM_H__66C50CA8_C861_41BC_A180_A2B3916CC96F__INCLUDED_)
