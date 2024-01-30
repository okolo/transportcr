// Stecker98IROSpectrum.h: interface for the CStecker98IROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STECKER98IROSPECTRUM_H__5A5EB53D_21C1_4DB0_894C_54103DB62A79__INCLUDED_)
#define AFX_STECKER98IROSPECTRUM_H__5A5EB53D_21C1_4DB0_894C_54103DB62A79__INCLUDED_

#include "Background.h"
#include "TableFunc.h"

/*
The spectrum is based on figure 1 of astro-ph/9808110
Z-dependence is not implemented
*/
class CStecker98IROSpectrum : public IBackgroundSpectrum
{
public:
	CStecker98IROSpectrum(double aZmax = 1000);
	virtual ~CStecker98IROSpectrum();
	virtual bool init();
	/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z);
	virtual double MaxZ() const { return m_Zmax; }
	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;
	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;
private:
	double m_Zmax;
	CVector m_X;
	CVector m_Y;
	CTableFunction* m_spectrum;
};

#endif // !defined(AFX_STECKER98IROSPECTRUM_H__5A5EB53D_21C1_4DB0_894C_54103DB62A79__INCLUDED_)
