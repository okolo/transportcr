#pragma once
#include "Coupling.h"
#include "EMCrosssec.h"

class CBackgroundTable;

namespace couplings{

class PP : public Coupling, protected CouplingParameters
{
	class Channel_gamma_e : public CouplingChannelT<PP>
	{
	public:
		Channel_gamma_e(PP* aCoupling, TParticle aSec) : CouplingChannelT<PP>(aCoupling, EPhoton, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_gamma_gamma : public CouplingChannelT<PP>
	{
	public:
		Channel_gamma_gamma(PP* aCoupling) : CouplingChannelT<PP>(aCoupling, EPhoton, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	PP(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	void DebugOutput(const char* aFolder) const;
protected:
	virtual void InitCoef(const CBackgroundTable& aBackgr) = 0;//PP part of former PropagCoef::InitICSandPP_V1()
	virtual void AdjustCoef() = 0;//PP part of former PropagCoef::LeeAdjustment()
	virtual void FinalCoefAdjustment();
	CVector	aph;
	CMatrix	ce;
	CVector	intPrdrPPlow;
	CVector	intPdrPPlow;
	EMCrosssec cs;
};

class PPold : public PP
{
public:
	PPold(void){};
protected:
	virtual void InitCoef(const CBackgroundTable& aBackgr);
	virtual void AdjustCoef();
};

class PPapr10 : public PP
{
public:
	PPapr10(void){};
protected:
	virtual void InitCoef(const CBackgroundTable& aBackgr);
	virtual void AdjustCoef();
};

}//namespace couplings{
