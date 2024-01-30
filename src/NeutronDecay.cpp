#include "NeutronDecay.h"
#include "Ranges.h"
#include <math.h>
#include "Units.h"
#include "ParticleData.h"

namespace couplings {

//neutron decay
const double NeutronDecay::NeutronDecayChannel::Neutron_lifetime=887.0; // (sec)
double NeutronDecay::NeutronDecayChannel::n_decay=0; // neutron decay probability in internal units
const double NeutronDecay::NeutronDecayChannel::delta_m_MeV=1.293318; //neutron & proton mass difference (MeV)
double NeutronDecay::NeutronDecayChannel::n_decay_k=0; // neutrino energy in CM frame (internal units)
double NeutronDecay::NeutronDecayChannel::n_decay_E=0; // electron energy in CM frame (internal units)

void NeutronDecay::NeutronDecayChannel::InitConstants()
{
	double Me = ParticleData::getParticleMass(EElectron);
// neutron decay related quantities
	n_decay=units.Tunit/Neutron_lifetime;
	double delta_m = delta_m_MeV/units.Eunit;
	n_decay_k=0.5*(delta_m-Me*Me/delta_m);
	n_decay_E=sqrt(n_decay_k*n_decay_k+Me*Me);
}

NeutronDecay::NeutronDecay(void)
{
	NeutronDecayChannel::InitConstants();

	AddChannel(new Channel_n_n(this));
	AddChannel(new Channel_n_p(this));
	AddChannel(new Channel_n_e_aen(this, EElectron));
	AddChannel(new Channel_n_e_aen(this, ENeutrinoAE));
}

void NeutronDecay::Channel_n_n::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CBinning& gamma = Ranges().midNucleonGamma();

	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, -n_decay/gamma[iPrim]);
	}
}

void NeutronDecay::Channel_n_p::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CBinning& gamma = Ranges().midNucleonGamma();

	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, n_decay/gamma[iPrim]);
	}
}

void NeutronDecay::Channel_n_e_aen::Coef(CMatrixAddOnlyView& aCoef) const
{
	
	const int nn = Ranges().nE();
	const CBinning& midGamma = Ranges().midNucleonGamma();
	const CBinning& midE = Ranges().midE();
	double s_const=BC().ss2-1.0/BC().ss2;
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		double gamma = midGamma[iPrim];
		double beta=sqrt(1.0-1.0/gamma/gamma);
		double one_minus_beta = 1./gamma/gamma;
		if (one_minus_beta>1e-5) {
			one_minus_beta = 1. - beta;
		}
		double decayRate = n_decay/gamma;
		double min_E, max_E;
		if(Secondary()==EElectron)
		{
			min_E=gamma*(n_decay_E-beta*n_decay_k);//electron min energy
			max_E=gamma*(n_decay_E+beta*n_decay_k);//electron max energy
		}
		else
		{//ENeutrinoAE
			min_E=gamma*n_decay_k*one_minus_beta;//neutrino min energy
			max_E=gamma*n_decay_k*(1.0+beta);//neutrino max energy
		}
		double diffSigmaNorm = 0.5*decayRate*s_const/(n_decay_k*beta*gamma);

		for(int iSec = 0; iSec<iPrim; iSec++)
		{
			double E1 = midE[iSec];
			if((E1>min_E)&&(E1<max_E))
				aCoef.Add(iSec, iPrim, diffSigmaNorm*E1);
		}
	}
}

}//namespace couplings {


