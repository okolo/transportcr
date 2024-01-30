// Stecker2005IROSpectrum.h: interface for the CStecker2005IROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(STECKER_2005_IRO__INCLUDED_)
#define STECKER_2005_IRO__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "TableBackground.h"
#include "TableFunc.h"

// IR/O spectrum from Stecker at al. astro-ph/0510449 
//(received from Matt Malkan <malkan@astro.UCLA.EDU> on the 6th of Jan 2006)

class CStecker2005IROSpectrum : public CTableWithHeaderBackgr  
{
public:
	CStecker2005IROSpectrum();
	virtual ~CStecker2005IROSpectrum();
// from CIROSpectrum
	double F(double E, double z);

	virtual double MaxZ() const;
	virtual double MinE(double aZmax) const;
	virtual double MaxE(double aZmax) const;

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

class Stecker2012IROSpectrum : public IBackgroundSpectrum{
public:
	Stecker2012IROSpectrum(bool aUpperLevel);
	virtual bool init();
	virtual double F(double E, double z);
	virtual double MaxZ() const;
	virtual double MaxE(double aZmax) const;
	virtual double MinE(double aZmax) const;
private:
	CStecker2005IROSpectrum f2005spec;
	MatrixBackground f2012specUV;
	double fMaxZ;
};

#endif // !defined(AFX_STECKER_2005_IRO__INCLUDED_)
