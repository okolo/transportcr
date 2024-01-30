#pragma once
#include "Coupling.h"

namespace couplings
{

class NeutronDecay : public Coupling
{
	class NeutronDecayChannel : public CouplingChannelT<NeutronDecay>
	{
	public:
		NeutronDecayChannel(NeutronDecay* aCoupling, TParticle aSec) : CouplingChannelT<NeutronDecay>(aCoupling, ENeutron, aSec){};
		static void InitConstants();
	protected:
		static const double Neutron_lifetime; // neutron lifetime in it's rest frame (sec)
		static double n_decay; // neutron decay probability it's rest frame (internal units)
		static const double delta_m_MeV; //neutron & proton mass difference (MeV)
		static double n_decay_k; // neutrino energy in CM frame (internal units)
		static double n_decay_E; // electron energy in CM frame (internal units)
	};

	class Channel_n_n : public NeutronDecayChannel
	{
	public:
		Channel_n_n(NeutronDecay* aCoupling) : NeutronDecayChannel(aCoupling, ENeutron){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
	class Channel_n_p : public NeutronDecayChannel
	{
	public:
		Channel_n_p(NeutronDecay* aCoupling) : NeutronDecayChannel(aCoupling, EProton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
	class Channel_n_e_aen : public NeutronDecayChannel
	{
	public:
		Channel_n_e_aen(NeutronDecay* aCoupling, TParticle aSec) : NeutronDecayChannel(aCoupling, aSec){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	NeutronDecay(void);
};

}//namespace couplings
