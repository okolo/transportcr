// SecondaryDifSigma.h: interface for the CSecondaryDifSigma class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SECONDARYDIFSIGMA_H__D633C233_51DB_11D5_B6C0_B76B625880B2__INCLUDED_)
#define AFX_SECONDARYDIFSIGMA_H__D633C233_51DB_11D5_B6C0_B76B625880B2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Sigma.h"
#include "DecaySpectrum.h"
#include "Addfunc.h"
#include "Vector.h"

class CSecondaryDifSigma : public CDiscretDifSigma  
{
public:
	double CalculateSigma_tot(int _s);
	virtual void Init();
	CSecondaryDifSigma(CDifSigma* _primarySigma, double _Smin, int _nn=Ranges().nE());
	virtual ~CSecondaryDifSigma();
	void AddSecondaryParticle(TParticle _particle);
	void AddPrimaryMode(TParticle _particle, double _rate = 1.);
	void AddDecayMode(DecaySpectrum* _mode, double _rate = 1.);

//creating links to avoid extra calculation, these functions must be called
//from PostInit() or after calling Init() function
	void AddSecondaryLink(TParticle _from, TParticle _to);//make secondary sigma of _from equal to one of _to
	void AddLink(TParticle _from, TParticle _to);//make sigma of _from equal to one of _to
	void MakeAntiparticlesLikeParticles();
	void DeleteSecondaryLink(TParticle _from);
	void MoveSecondaryLink(TParticle _from, TParticle _to);


//implemented from CDiscretDifSigma
	inline double Sigma(int _s,int _ratio,TParticle _particle) const;//sigma multiplied by E (initial particle energy)
protected:
	virtual void PostInit(){};//called at the end of Init
	class TDecayMode
	{
	public:
		DecaySpectrum* iSpectrum;
		double iRate;
		TDecayMode(DecaySpectrum* aSpectrum=NULL, double aRate=0):
			iSpectrum(aSpectrum),iRate(aRate){};
			~TDecayMode(){/*delete iSpectrum;*/};
	};
	CDifSigma* m_difSigma;
	CVector m_directRates;
	int m_maxNoOfDecayModes;//m_secondaryModes array size
	int m_noOfDecayModes;//actual number of nonzero elements in m_secondaryModes array
	bool isReady;
	int m_NoOfSigmasMem;//the actual m_sigmasMem array size (<=EEndParticle)
	CPointerArray<CMatrix> m_sigmas;//links to calculated sigmas (m_sigmasMem) with order corresponding to TParticle enum
	CAutoDeletePtrArray<CMatrix> m_sigmasMem;//actual array of matrixes in undefined order
	CMatrix m_primarySigma;//primary chanel sigma
	TDecayMode* m_secondaryModes;
private:
	bool isPreLinkInitDone;//flag indicating if all the initialization exept
	//establishing links has been done
};

inline double CSecondaryDifSigma::Sigma(int _s,int _ratio,TParticle _particle) const
{
	ASSERT(isReady);
	ASSERT(_s<m_nn&&_ratio<Ranges().nE());
	double result = m_directRates[_particle]*m_primarySigma[_s][_ratio];
	if((m_sigmas[_particle]!=NULL)&&(_ratio<Ranges().nE()-1))
		result+=m_sigmas(_particle)[_s][_ratio];
	ASSERT_VALID_NO(result);
	return result;
}


#endif // !defined(AFX_SECONDARYDIFSIGMA_H__D633C233_51DB_11D5_B6C0_B76B625880B2__INCLUDED_)
