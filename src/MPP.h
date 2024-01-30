#pragma once
#include "Coupling.h"
#include "MuPP.h"
#include "Vector.h"
#include "const.h"
#include "EMCrosssec.h"

namespace couplings{

//Muon-antimuon pair production by photons
class MPP : public Coupling, CouplingParameters, EMConstants
{
	class Channel_gamma_e : public CouplingChannelT<MPP>
	{
	public:
		Channel_gamma_e(MPP* aCoupling, TParticle aSec) : CouplingChannelT<MPP>(aCoupling, EPhoton, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_gamma_gamma : public CouplingChannelT<MPP>
	{
	public:
		Channel_gamma_gamma(MPP* aCoupling) : CouplingChannelT<MPP>(aCoupling, EPhoton, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_gamma_neutrino_e : public CouplingChannelT<MPP>
	{
	public:
		Channel_gamma_neutrino_e(MPP* aCoupling, TParticle aSec) : CouplingChannelT<MPP>(aCoupling, EPhoton, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_gamma_neutrino_mu : public CouplingChannelT<MPP>
	{
	public:
		Channel_gamma_neutrino_mu(MPP* aCoupling, TParticle aSec) : CouplingChannelT<MPP>(aCoupling, EPhoton, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	MPP(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
private:
	double RmuPP(double E,double b);

	CMatrix g_e;
	CMatrix g_mn;
	CMatrix g_en;
	CVector a_ph;

	EMCrosssec cs;
};

}//namespace couplings{
