
#include <fstream>
#include "Addfunc.h"
#include <math.h>
#include "Ranges.h"
#include "Nucleus.h"
#include "TimeZ.h"
#include "Units.h"

extern const double alfa_1;

#include "Ppp.h"
#include "Ranges.h"
#include "Medium.h"

namespace couplings{

//energy loss rate
double PPP::ProtonFi(double _k) const
{// all the local calculations in system of units Me=1
	const double c1=0.8048;
	const double c2=0.1459;
	const double c3=1.137e-3;
	const double c4=-3.879e-6;
	const double d0=-86.07;
	const double d1=50.96;
	const double d2=-14.45;
	const double d3=2.6666667;
	const double f1=2.910;
	const double f2=78.35;
	const double f3=1837;
	const double Pi_12=0.26179939;
	double result;
	if(_k<=2.0)
		return 0.0;
	if(_k<25.0)
	{
		double k=_k-2;
		result=Pi_12*k*k*k*k/(1.0+k*(c1+k*(c2+k*(c3+k*c4))));
	}
	else
	{
		double lnk=log(_k);
		result=_k*(d0+lnk*(d1+lnk*(d2+lnk*d3)));
		result/=(1.0-(f1+(f2+f3/_k)/_k)/_k);
	};
	return result/_k/_k;//*1e5;//1e5 artifitial coef. for testing energy conserv.
}

/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  /*  Pair Production by Protons (PPP) angle averaged differential
  cross section (per e energy) multiplied by E  */
double PPP::P_PPP(double E,double b,double nE) const
/* E-cosmic ray proton energy
   b-background photon energy
   nE-the energy of one of the particles of the produced pair */
{
    ASSERT(cPPP1>0.);//constants must be initialized
	double eMax,eMin,s,result;
	double p = 0.;
	double beta = 0.;
	if (E<=10*Mp){// treat E parameter as momentum
        p = E;
        E = sqrt(Mp2 + p*p);
	}
	else{
        p = sqrt(E*E-Mp2);
	}
    beta = p/E;
	s = E*b;
	if(s*(1.+beta)<thresholdPPP)
		return 0.;

	const double delta = 1.75;//  7/4,  P ~ nE^{-delta}
	double d1 = 2. - delta;
	double r = nE/E;

	//obtaining eMax & eMin formulas were taken from Lee article
	// should be proved !!!


	double s0 = b*(E + p);
	double x = s0*(s0-2.*Me*(Me+lambdaPPP))-cPPP1;
	if(x<0.)
		return 0.;
	double dif = (b+p)*sqrt(x);
	double middle = (E+b)*(s0-Me*lambdaPPP);
	double mult = 1./(2.*s0+Mp2);
	if(middle<=dif)
	{//TODO (see pragma below)
//#pragma message (__FILE__ "(" STRING(__LINE__) "): todo: PPP make diff. cross section independent on Emin (investigate accuracy problems)")
		//ASSERT(middle==dif);
		eMin = Ranges().Emin();
	}
	else
		eMin = (middle - dif)*mult;
	eMax = (middle + dif)*mult;

	if((nE>eMax)||(nE<eMin))
		return(0.);

	double rMin = eMin/E;
	double rMax = eMax/E;
	double c = 0.5*d1/(pow(rMax,d1)-pow(rMin,d1));//0.5 multiplier comes from dividing energy between e+ and e-
	// obtained from energy conservation in the same way
	//as for TPP

	result = PPP_const*2./Mp*ProtonFi(2.*b*E/Mp/Me);
	//double oldResult =  result*(0.125*E/(pow(eMax,0.25)-pow(eMin,0.25)))*pow(nE,-1.75);
	result *= (pow(r,-delta)*c /* 1.283834 */);
	return result;
}

bool PPP::NO_PPP_PRODUCTS = false;

PPP::PPP(void):
	Mp(ParticleData::getParticleMass(EProton)),
	Mp2(Mp*Mp),
	cPPP1(0.5*Me2*(Mp2-Me2)),
	lambdaPPP(sqrt(0.5*(Mp2+Me2))),
	thresholdPPP(2*Me*(Mp+Me)),
	PPP_const(alpha*alpha*alpha/Me)
{
	READ_BOOL_SETTING(NO_PPP_PRODUCTS);
	AddChannel(new Channel_A_A(this, EProton));
	if(!NO_PPP_PRODUCTS)
	{
		AddChannel(new Channel_A_e(this, EProton, EElectron));
		AddChannel(new Channel_A_e(this, EProton, JOIN_ELECTRON_POSITRON ? EElectron : EPositron));	
	}

	FOR_ALL_NUCLEI_INVOLVED(particle)
	{
		AddChannel(new Channel_A_A(this, particle));
		if(!NO_PPP_PRODUCTS)
		{
			AddChannel(new Channel_A_e(this, particle, EElectron));
			AddChannel(new Channel_A_e(this, particle, EPositron));
		}
	}
}

void PPP::SetBackgrounds(const Medium& aPropagCoef)
{
	Coupling::SetBackgrounds(aPropagCoef);

	const int nn = Ranges().nE();
	if(!appp.length())
	{
		appp.create(nn);
		pe_PPP.create(nn);
	}
	else
	{
		appp.reset();
		pe_PPP.reset();
	}
	const SmartPtr<CBackgroundTable>& ebl = aPropagCoef.background();
	const int mm = ebl->nK();
	const double s_const=BC().ss2-1.0/BC().ss2;
	const double p_mass = ParticleData::getParticleMass(EProton);
	for(int k=0; k<mm; k++)
	{
		double b = ebl->midBackgroundE()[k];
		double Fb=ebl->Fmid(k);
		for(int i=0; i<nn; i++)
		{
			appp[i]+=Fb*PPP_const*2.0/p_mass*
				ProtonFi(2*b*Ranges().E()[i]/p_mass/Me)/s_const;//appp[i]=1/E*dE/dt/(ss1-1)
		};
	};
	if(NO_PPP_PRODUCTS)
		return;

	const CBinning& midE = Ranges().midE();
	const CBinning& leftE = Ranges().E();
	for(int j=0; j<nn; j++)
	{
		double E = midE[j];

		for(int k=0; k<mm; k++)
		{
			double b = ebl->midBackgroundE()[k];
			double sB_const=s_const*ebl->Fmid(k);

			for(int i=0; i<j; i++)
			{
				double E1 = midE[i];
				double factor=sB_const*E1/E;
				double add = P_PPP(leftE[j],b, E1)*factor;
                ASSERT_VALID_NO(add);
				pe_PPP[i][j] +=	add;
			};
		};
	};
}

void PPP::DebugOutput(const char* aFolder) const
{
// print energy loss rate -dE/dt for protons
// real E loss rate is given by integral ProtonFi n(b)db
// effective rate is calculated based on proton propagation coefficients
// discrepancy between these two rates may mean that something is wrong with current CEL simulation method
// at the time of writing this text work on PPP and TPP accuracy enhancement was in progress in branch tpp_ppp_check
	const CBackgroundTable* ebl = Background()->background();
	const int nn = Ranges().nE();
	const int mm = ebl->nK();
	string path = aFolder;
	path = path + "PPP";
	Mkdir(path);
	path = path + DIR_DELIMITER_STR + ToString(redshift.z());
	ofstream file(path.c_str());
	const CBinning& E=Ranges().E();//using Ranges().E and not Ranges().midE for CEL approximation
	const CBinning& b=ebl->midBackgroundE();
	const double p_mass = ParticleData::getParticleMass(EProton);

	double dt_dz = CTimeZ::diffT(redshift.z());
	double s_const=BC().ss2-1.0/BC().ss2;

	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		double eLoss = 0.;
		for(int k=0; k<mm; k++)
		{
			double Fb=ebl->Fmid(k);
			eLoss += ProtonFi(2*b[k]*E[iPrim]/p_mass/Me)*Fb;
		}
		eLoss *= (PPP_const*2.0/p_mass*E[iPrim]*units.Eunit*1e6/units.Lunit*units.Mpc_cm);//[eV/Mpc]
		double eLossZ = eLoss*units.Lunit/units.Mpc_cm *dt_dz;// dE/dz [eV]
		double eLossEffective = E[iPrim]*appp[iPrim]*s_const*units.Eunit*1e6/units.Lunit*units.Mpc_cm;//[eV/Mpc]
		if(eLossEffective>0. || eLoss>0.)
			file << (E[iPrim]*units.Eunit*1e6) << "\t" << eLossZ << "\t" << eLoss << "\t" << eLossEffective << "\n";
	}
	file.close();
}

void PPP::GetMultPPP(double& multP, double& multR, TParticle aPrim)
{
	double Z = CNucleus::getZ(aPrim);
	double A = CNucleus::getA(aPrim);
	multP = Z*Z;
	multR = multP/((double)A);
}

void PPP::Channel_A_A::Coef(CMatrixAddOnlyView& aCoef) const
{
	double multR, multP;
	PPP::GetMultPPP(multP, multR, Primary());
	//const PropagCoef& backgr = *Background();
	const int nn = Ranges().nE();
	int iPrim;
	for(iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, - multR * fCoupling.appp[iPrim]);
		if(iPrim<nn-1)
			aCoef.Add(iPrim, iPrim+1, multR * fCoupling.appp[iPrim+1]);
	}
}

void PPP::Channel_A_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	double multR, multP;
	PPP::GetMultPPP(multP, multR, Primary());
	//const PropagCoef& backgr = *Background();
	const int nn = Ranges().nE();
	int iPrim, iSec;
	for(iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		for(iSec = 0; iSec<=iPrim; iSec++){
            ASSERT_VALID_NO(fCoupling.pe_PPP[iSec][iPrim]);
            aCoef.Add(iSec, iPrim, multP*fCoupling.pe_PPP[iSec][iPrim]);
		}
	}
}

}//end of namespace couplings

