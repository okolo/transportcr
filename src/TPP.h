#pragma once
#include "Coupling.h"
#include "const.h"

namespace couplings
{
class TPP : public Coupling, CouplingParameters, EMConstants
{
	class Channel_e_e : public CouplingChannelT<TPP>
	{
	public:
		Channel_e_e(TPP* aCoupling, TParticle aPrim) : CouplingChannelT<TPP>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
	class Channel_e_ae : public CouplingChannelT<TPP>
	{
	public:
		Channel_e_ae(TPP* aCoupling, TParticle aPrim, TParticle aSec) : CouplingChannelT<TPP>(aCoupling, aPrim, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	TPP(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	virtual void DebugOutput(const char* aFolder) const;
private:
	double RTPP(double E,double b);
	double F_TPP(double E,double b);
	double PTPP(double E,double b,double nE);

	CVector f_a;
	CMatrix f_b;
};
}
