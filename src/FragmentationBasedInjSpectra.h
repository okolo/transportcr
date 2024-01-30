#ifndef FragmentationBasedInjSpectra_H_ALREADY_INCLUDED
#define FragmentationBasedInjSpectra_H_ALREADY_INCLUDED


#include "InjectionSpectra.h"
#include "Vector.h"
#include "TableFunc.h"

class CFragmentationBasedInjSpectra :
	public CInjectionSpectra
{
public:
	CFragmentationBasedInjSpectra();
	virtual ~CFragmentationBasedInjSpectra(void);
	virtual double Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ);
	virtual void init();


private:
	CAutoDeletePtrArray<CTableFunction>	m_fragSpectra;
	double								m_graphPower;//power base of graph data (default is 3, that is using x^3*F(x) graph as source)
	const char*							m_SourceDir;
};

#endif //#ifndef FragmentationBasedInjSpectra_H_ALREADY_INCLUDED
