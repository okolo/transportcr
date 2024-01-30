// UserInjSpectra.h: interface for the CTableInjSpectra class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_USERINJSPECTRA_H__5F5DD021_C191_11D5_885E_444553540001__INCLUDED_)
#define AFX_USERINJSPECTRA_H__5F5DD021_C191_11D5_885E_444553540001__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "InjectionSpectra.h"
#include "Concentrations.h"
#include "Vector.h"
#include "ParticleData.h"
#include "Parameters.h"

class CTableInjSpectra : public CInjectionSpectra
{
public:
	CTableInjSpectra();
	//CTableInjSpectra(const Concentrations& conc, double aDeltaT_Mpc);
	virtual double Q(TParticle aParticle, int aBinE, double aE, double aZ);
	virtual void init();
protected:
	SafePtr<Concentrations>		m_spectra;
	std::string					m_DataFileName;
	double						m_deltaT;//internal units
};

#endif // !defined(AFX_USERINJSPECTRA_H__5F5DD021_C191_11D5_885E_444553540001__INCLUDED_)
