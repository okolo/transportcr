#include "MPP.h"
#include "Ranges.h"
#include "Medium.h"
#include "Units.h"

extern double RPP(double _E,double _b);
namespace couplings{

double MPP::RmuPP(double E,double b)
{
	return cs.RPP(E,b);
}

MPP::MPP(void):
cs(Mmu)
{
	AddChannel(new Channel_gamma_gamma(this));
	AddChannel(new Channel_gamma_e(this, EElectron));
	AddChannel(new Channel_gamma_e(this, JOIN_ELECTRON_POSITRON ? EElectron : EPositron));
	AddChannel(new Channel_gamma_neutrino_e(this, ENeutrinoE));
	AddChannel(new Channel_gamma_neutrino_e(this, ENeutrinoAE));
	AddChannel(new Channel_gamma_neutrino_mu(this, ENeutrinoM));
	AddChannel(new Channel_gamma_neutrino_mu(this, ENeutrinoAM));
}

void MPP::Channel_gamma_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		aCoef.Add(iPrim, iPrim, -fCoupling.a_ph[iPrim]);
}

void MPP::Channel_gamma_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim, iSec;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, fCoupling.g_e[iSec][iPrim]);
	}
}

void MPP::Channel_gamma_neutrino_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim, iSec;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, fCoupling.g_en[iSec][iPrim]);
	}
}

void MPP::Channel_gamma_neutrino_mu::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim, iSec;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, fCoupling.g_mn[iSec][iPrim]);
	}
}

void MPP::SetBackgrounds(const Medium& aPropagCoef)
{
	Coupling::SetBackgrounds(aPropagCoef);
	SmartPtr<CBackgroundTable> background = aPropagCoef.background();
	CMuPairProduction MuPairProduction(background->Kmin(), background->nK());
	const int nn = Ranges().nE();
	if(!a_ph.length())
	{
		a_ph.create(nn);
		g_e.create(nn);
		g_mn.create(nn);
		g_en.create(nn);
	}
	else
	{
		a_ph.reset();
		g_e.reset();
		g_mn.reset();
		g_en.reset();
	}
	double Emin = Ranges().Emin();
	double ss2 = BC().ss2;
	double ss1 = BC().ss1;
	double s_const=BC().ss2-1.0/BC().ss2;
	int mm = background->nK();
	int i,j,k;
	double E;
	for(j=0, E=Emin*ss2;j<nn; j++, E*=ss1)
	{
		for(k=0; k<mm; k++)
		{	
			double Fb=background->Fmid(k);
			a_ph[j] += Fb*RmuPP(E, background->midBackgroundE()[k]);
			double sB_const=s_const*Fb;
			double E1;
			for(i=0,E1=Emin*ss2;i<j;i++,E1*=ss1)
			{
				double factor=sB_const*E1/E;
				g_en[i][j]+=(MuPairProduction.GetPen(i, j, k)*factor);
				g_mn[i][j]+=(MuPairProduction.GetPmn(i, j, k)*factor);
				g_e[i][j]+=(MuPairProduction.GetPe(i, j, k)*factor);
			};
		};
	};
}

}//namespace couplings{

