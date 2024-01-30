#pragma once

#include "TableFunc.h"
#include "Vector.h"
#include "Background.h"
#include "Coupling.h"
#include "const.h"

namespace couplings{

class Dpp : public Coupling, IFunction, CouplingParameters, EMConstants
{
	class Channel_gamma_e : public CouplingChannelT<Dpp>
	{
	public:
		Channel_gamma_e(Dpp* aCoupling, TParticle aSec) : CouplingChannelT<Dpp>(aCoupling, EPhoton, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_gamma_gamma : public CouplingChannelT<Dpp>
	{
	public:
		Channel_gamma_gamma(Dpp* aCoupling) : CouplingChannelT<Dpp>(aCoupling, EPhoton, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	Dpp();
	~Dpp(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	inline double bDPP(int aElectronE, int aPhotonE) const
	{
		double r =  m_midE[aElectronE]/m_midE[aPhotonE];
		//return m_sF * r * m_fi->f(r) * m_SigmaBackgroundIntegral[aPhotonE];  //m_Backgr.integralDpp(m_threshold/m_midE[aPhotonE]);
		return m_sF * r *analyticalFi(r)* m_SigmaBackgroundIntegral[aPhotonE];//analytical approximation appears to be accurate enough
	}

	inline double aDPP(int aPhotonE) const
	{
		return 2.*m_SigmaBackgroundIntegral[aPhotonE];
	}

	inline double analyticalFi(double aR) const
	{
		if(aR<0 || aR>1)
			return 0.;
		double dr = aR-0.5;
		return (1.66666666666667 + (4.*dr*dr));//based on the parametrization of the table function
	}

	
private:
	void SetBackground(const CBackgroundIntegral& aBackgr);
	virtual double f(double _x) const;//used for background integral calculation
	double m_fParamE;
	void InitSigmaTot();
	double			m_sF;
	double			m_threshold;
	CLinearFunc*	m_fi;
	CLinearFunc*	m_SigmaTot;
	CVector			m_SigmaBackgroundIntegral;
	CVector			m_x;
	CVector			m_y;
	const CBackgroundIntegral* m_Backgr;
	const CBinning&	m_midE;
};

}//namespace couplings
