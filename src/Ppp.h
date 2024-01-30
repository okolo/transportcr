#pragma once
#include "Coupling.h"
#include "const.h"

namespace couplings{

class PPP : public Coupling, CouplingParameters, EMConstants
{
	class Channel_A_e : public CouplingChannelT<PPP>
	{
	public:
		Channel_A_e(PPP* aCoupling, TParticle aPrim, TParticle aSec) : CouplingChannelT<PPP>(aCoupling, aPrim, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_A_A : public CouplingChannelT<PPP>
	{
	public:
		Channel_A_A(PPP* aCoupling, TParticle aPrim) : CouplingChannelT<PPP>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	PPP(void);
	static void GetMultPPP(double& multP, double& multR, TParticle aPrim);
	virtual void DebugOutput(const char* aFolder) const;
	virtual void SetBackgrounds(const Medium& aPropagCoef);
private:
	double ProtonFi(double _k) const;
	double P_PPP(double E,double b,double nE) const;

	CVector appp;
	CMatrix pe_PPP;
	const double Mp;
	const double Mp2;
	const double cPPP1;
	const double lambdaPPP;
	const double thresholdPPP;
	const double PPP_const;
	static bool NO_PPP_PRODUCTS;
};

}//namespace couplings{
