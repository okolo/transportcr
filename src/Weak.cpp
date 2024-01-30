#include "Weak.h"
#include "Ranges.h"
#include "Medium.h"
#include "const.h"
#include "Units.h"

typedef CSecondaryDifSigma* CPSecondaryDifSigma;
typedef CNeutrinoSNonHadronChanel* CPNeutrinoSNonHadronChanel;
namespace couplings{



/*
bool NO_NEUTRINO_t;
bool NO_NEUTRINO_s_HADRON);
bool NO_NEUTRINO_s_NON_HADRON);
NO_NEUTRINO_s = NO_NEUTRINO_s_HADRON&&NO_NEUTRINO_s_NON_HADRON;

bool Z_CHANEL_YOSHIDA;
bool JOIN_ELECTRON_POSITRON;
double NEUTRINO_S_CHANNEL_MULTIPLIER;
double NEUTRINO_T_CHANNEL_MULTIPLIER;
int RESONANCE_RESOLUTION_ACCURACY;*/



WeakCoupling::~WeakCoupling(void)
{
	Close();
}

void WeakCoupling::Init()
{
	const double Emin = Ranges().Emin();
	const int nn = Ranges().nE();
	if(m_neutrinoChanals)
		return;

	FWeak::Init();
	CZHadronChanel::InitLEPData();

	m_neutrinoChanals = new CPSecondaryDifSigma[end_tE];
	int i;
	for(i=0;i<end_tE;i++)
		m_neutrinoChanals[i]=0;

	m_sNonHadronChanels = new CPNeutrinoSNonHadronChanel[end_sE];
	m_MassiveNeutrinoSHadronChanels = new CMassiveNeutrinoSHadronChanel[end_sE];
	for(i=0;i<end_sE;i++)
	{
		m_sNonHadronChanels[i]=0;
	}
	
	if(!NO_NEUTRINO_t)
	{
		m_neutrinoChanals[enaenE]= new CEnAenChanel(2.*m_en*Emin*BC().ss2);
		m_neutrinoChanals[enamnbE]= new CEnAmnChanel(2.*m_mn*Emin*BC().ss2);
		m_neutrinoChanals[enbamnE]= new CEnAmnChanel(2.*m_en*Emin*BC().ss2);
		m_neutrinoChanals[enatnbE]= new CEnAtnChanel(2.*m_tn*Emin*BC().ss2);
		m_neutrinoChanals[mnamnE]= new CMnAmnChanel(2.*m_mn*Emin*BC().ss2);
		m_neutrinoChanals[mnatnbE]= new CMnAmnChanel(2.*m_tn*Emin*BC().ss2);
		for(i=0;i<end_tE;i++)
			m_neutrinoChanals[i]->Init();
		

		m_tChanelSigma = (CNeutrinoTChanelSigma*)(m_neutrinoChanals[0]);
//#ifdef _DEBUG
		m_tChanelSigma->Print(plt_local_dir + "t_tot",2.*m_en*Emin,nn);
		m_tChanelSigma->printSigmaTot(plt_local_dir + "t_dif",2.*m_en*Emin,nn,2*nn);
		m_tChanelSigma->printSigmaTotE(plt_local_dir + "t_difE",2.*m_en*Emin,nn,2*nn);
//#endif
	}
	
	if(!NO_NEUTRINO_s_NON_HADRON)
	{
		m_sNonHadronChanels[s_eE]= new CNeutrinoSNonHadronChanel(2.*m_en*Emin*BC().ss2);
		m_sNonHadronChanels[s_mE]= new CNeutrinoSNonHadronChanel(2.*m_mn*Emin*BC().ss2);
		m_sNonHadronChanels[s_tE]= new CNeutrinoSNonHadronChanel(2.*m_tn*Emin*BC().ss2);
		for(i=0;i<end_sE;i++)
			m_sNonHadronChanels[i]->Init();
	}
	
	if(!NO_NEUTRINO_s)// even if NO_NEUTRINO_s_HADRON = true we still need
	//to calculate m_MassiveNeutrinoSHadronChanels (see K_sMNChanel)
	{
		m_MassiveNeutrinoSHadronChanels = new CMassiveNeutrinoSHadronChanel[end_sE];
		m_MassiveNeutrinoSHadronChanels[s_eE].Init(m_en);
		m_MassiveNeutrinoSHadronChanels[s_mE].Init(m_mn);
		m_MassiveNeutrinoSHadronChanels[s_tE].Init(m_tn);
		m_MassiveNeutrinoSHadronChanels[s_eE].Print(plt_local_dir + "s_tot",2.*m_en*Emin,nn);
	}	
}

//Former PropagCoef::CloseNeutrinoChanals()
void WeakCoupling::Close()
{
	if(!m_neutrinoChanals)
		return;
	int i;
	for(i=0;i<end_tE;i++)
		delete m_neutrinoChanals[i];
	delete[] m_neutrinoChanals;
	for(i=0;i<end_sE;i++)
		delete m_sNonHadronChanels[i];
	delete[] m_sNonHadronChanels;

	delete[] m_MassiveNeutrinoSHadronChanels;
	
	m_neutrinoChanals = 0;
	m_tChanelSigma = 0;
	m_sNonHadronChanels = 0;
	m_MassiveNeutrinoSHadronChanels = 0;

	CZHadronChanel::CloseLEPData();
}

//Former PropagCoef::WriteNeutrinoCoef()
void WeakCoupling::WriteCoef()
{
	const double Emin = Ranges().Emin();
	const int nn = Ranges().nE();
	if(NO_NEUTRINO_s&&NO_NEUTRINO_t)
		return;
	double E = Emin*BC().ss2*units.Eunit*1e6;// eV
	int i = 0;
	double mult = units.Lunit/units.Mpc_cm /1000;//Gpc
	double neutrinoConc = fBackground->neutrinoConc;
	double sigmaRes = m_MassiveNeutrinoSHadronChanels[s_eE].SigmaTotalResonanceValue();
	double resonanceL = mult/neutrinoConc/sigmaRes;
	sigmaRes*=units.sigmaUnit;//sigma in cm^2

	FILE* file = Fopen(plt_local_dir + "neutrinoIntL","wt");
	fprintf(file,"#background neutino concentration (all 6 types) %lg cm^-3\n",neutrinoConc/units.Vunit*6.0);//6 neutino types
	fprintf(file,"#Z-resonance: sigma = %lg cm^2, L_int = %lg Gpc\n",sigmaRes,resonanceL);
	fprintf(file,"#electron neutrino interaction lengths (Gpc)\n"
	"# E,eV\t\t s-channel\t\t t-channel\n");
	for(;i<nn;i++,E*=BC().ss1)
	{//todo: port coefficient output procedure
		double tL = 0;// mult/Concentrations::K_tMNChanel(0,*this,i);
		double sL = 0;// mult/Concentrations::K_sMNChanel(0,*this,i,s_eE);
		fprintf(file,"%lg\t%lg\t%lg\n",E,sL,tL);
	}
	fclose(file);
}



WeakCoupling::TNeutrinoSChanels WeakCoupling::GetSChanelIndex(TParticle aPrimary)
{
	switch(aPrimary)
	{
	case ENeutrinoAE:
	case ENeutrinoE:
		return s_eE;
	case ENeutrinoAM:
	case ENeutrinoM:
		return s_mE;
	case ENeutrinoAT:
	case ENeutrinoT:
		return s_tE;
	default:
	  ASSERT(0);
	  throw "invalid argument passed to WeakCoupling::GetSChanelIndex";
	}
	return end_sE;
}

WeakCoupling::Channel_HadronS::Channel_HadronS(WeakCoupling* aCoupling, TParticle aPrim, TParticle aSec):
CouplingChannelT<WeakCoupling>(aCoupling, aPrim, aSec),
fChannel(aCoupling->m_MassiveNeutrinoSHadronChanels[WeakCoupling::GetSChanelIndex(aPrim)])
{
	switch(aSec)
	{
	case EElectron:
	case EPositron:
		fSecondaryGroup = CLEPData::electronsE;
		fMult = 1./2.;//two secondaries
		break;
	case EPhoton:
		fSecondaryGroup = CLEPData::photonsE;
		fMult = 1.;//one secondary
		break;
	case ENeutron:
		fSecondaryGroup = CLEPData::neutronsE;
		fMult = 1.;//one secondary
		break;
	case EProton:
		fSecondaryGroup = CLEPData::protonsE;
		fMult = 1.;//one secondary
		break;
	case ENeutrinoAE:
	case ENeutrinoAM:
	case ENeutrinoAT:
	case ENeutrinoE:
	case ENeutrinoM:
	case ENeutrinoT:
		fSecondaryGroup = CLEPData::neutrinosE;
		fMult = 1./6.;//6 secondaries
		break;
	default:
		ASSERT(0);
	}
}

void WeakCoupling::Channel_HadronS::Coef(CMatrixAddOnlyView& aCoef) const
{
	const Medium& backgr = *Background();
	double neutrinoConc =  backgr.neutrinoConc*backgr.neutrinoClusteringModifier;
	const int nn = Ranges().nE();
	double constMult = (BC().ss2-BC().ss_2)*neutrinoConc*fMult;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		double mult = constMult * Ranges().midE()[iPrim];
		for(int iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, mult*fChannel.P(iPrim, iSec, fSecondaryGroup));
	}
}

WeakCoupling::Channel_NonHadronS::Channel_NonHadronS(WeakCoupling* aCoupling, TParticle aPrim, TParticle aSec):
CouplingChannelT<WeakCoupling>(aCoupling, aPrim, aSec),
fChannel(*aCoupling->m_sNonHadronChanels[WeakCoupling::GetSChanelIndex(aPrim)])
{
}
void WeakCoupling::Channel_NonHadronS::Coef(CMatrixAddOnlyView& aCoef) const
{
	const Medium& backgr = *Background();
	double neutrinoConc =  backgr.neutrinoConc*backgr.neutrinoClusteringModifier;
	const int nn = Ranges().nE();
	const CBinning& E = Ranges().midE();
	double constMult = (BC().ss2-BC().ss_2)*neutrinoConc;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		double mult = constMult / E[iPrim];
		for(int iSec = 0; iSec<=iPrim; iSec++)
		{
			aCoef.Add(iSec, iPrim, mult*E[iSec]*fChannel.Sigma(iPrim, nn-1-(iPrim-iSec), Secondary()));
		}
	}
}

void WeakCoupling::Channel_AbsorptionS::Coef(CMatrixAddOnlyView& aCoef) const
{
	const Medium& backgr = *Background();
	double mult = fMult* backgr.neutrinoConc*backgr.neutrinoClusteringModifier;
	const int nn = Ranges().nE();

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, - mult*fChannel.R(iPrim));
	}
}

WeakCoupling::Channel_AbsorptionS::Channel_AbsorptionS(WeakCoupling* aCoupling, TParticle aPrim):
CouplingChannelT<WeakCoupling>(aCoupling, aPrim, aPrim),
fChannel(aCoupling->m_MassiveNeutrinoSHadronChanels[WeakCoupling::GetSChanelIndex(aPrim)])
{
	fMult =  1.;
	if(NO_NEUTRINO_s_HADRON)
		fMult = 1. - FWeak::gOut_q;
	else if(NO_NEUTRINO_s_NON_HADRON)
		fMult = FWeak::gOut_q;
}

void WeakCoupling::Channel_AbsorptionT::Coef(CMatrixAddOnlyView& aCoef) const
{
	const Medium& backgr = *Background();
	double neutrinoConc = backgr.neutrinoConc*backgr.neutrinoClusteringModifier;
	const int nn = Ranges().nE();
	const CNeutrinoTChanelSigma& sigma = *fCoupling.m_tChanelSigma;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, - neutrinoConc*sigma.SigmaTotMassive(iPrim));
	}
}

WeakCoupling::Channel_T::Channel_T(WeakCoupling* aCoupling, const WeakCoupling::TChannelMapEntry& aEntry, TParticle aSec):
CouplingChannelT<WeakCoupling>(aCoupling, aEntry.Primary, aSec),
fChannel(*aCoupling->m_neutrinoChanals[aEntry.Channel])
{
	fEffSecondary = aEntry.UseAntiparticleCoef ? ParticleData::GetAntiparticle(aSec) : aSec;
}

WeakCoupling::WeakCoupling(void)
{
	if(JOIN_ELECTRON_POSITRON)
	{
		throw "JOIN_ELECTRON_POSITRON mode is not supported yet by WeakCoupling coupling";
	}
	const CParticleList& pf = *CParticleList::Instance();
	Init();
	
	if(!NO_NEUTRINO_s)
	{
		for(int primary = EStartNeutrino; primary < EEndNeutrino; primary++)
		{
			if(!pf.IsEnabled((TParticle)primary)) continue;
			AddChannel(new Channel_AbsorptionS(this, (TParticle)primary));
			if(!NO_NEUTRINO_s_HADRON)
			{
				for(int secondary = EStartEM; secondary <= EProton; secondary++)
					if(pf.IsEnabled((TParticle)secondary))
						AddChannel(new Channel_HadronS(this, (TParticle)primary, (TParticle)secondary));
			}
			if(!NO_NEUTRINO_s_NON_HADRON)
			{
				const TParticle secondaries[] = {EElectron,EPositron,ENeutrinoE,ENeutrinoM,
					ENeutrinoT,ENeutrinoAE,ENeutrinoAM,ENeutrinoAT,		EEndParticle};
				for(int i=0; secondaries[i]!=EEndParticle; i++)
					if(pf.IsEnabled(secondaries[i]))
						AddChannel(new Channel_NonHadronS(this, (TParticle)primary, secondaries[i]));
			}
		}
	}
	
	if(!NO_NEUTRINO_t)
	{
		const TParticle secondaries[] = {EElectron,EPositron,ENeutrinoE,ENeutrinoM,
				ENeutrinoT,ENeutrinoAE,ENeutrinoAM,ENeutrinoAT,		EEndParticle};
		for(int ic=0; fTChannelMap[ic].Primary != EEndParticle; ic++)
		{
			if(pf.IsEnabled(fTChannelMap[ic].Primary))
				for(int i=0; secondaries[i]!=EEndParticle; i++)
					if(pf.IsEnabled(secondaries[i]))
						AddChannel(new Channel_T(this, fTChannelMap[ic], secondaries[i]));
		}
		for(int particle=EStartNeutrino; particle < EEndNeutrino; particle++)
			if(pf.IsEnabled((TParticle)particle))
				AddChannel(new Channel_AbsorptionT(this, (TParticle)particle));
	}
}

const WeakCoupling::TChannelMapEntry WeakCoupling::fTChannelMap[] =
{
	{ENeutrinoE, enaenE, false},
	{ENeutrinoE, enamnbE, false},
	{ENeutrinoE, enatnbE, false},
	{ENeutrinoAE, enaenE, false},
	{ENeutrinoAE, enamnbE, true},
	{ENeutrinoAE, enatnbE, true},

	{ENeutrinoM, mnamnE, false},
	{ENeutrinoM, mnatnbE, false},
	{ENeutrinoM, enbamnE, true},
	{ENeutrinoAM, mnamnE, false},
	{ENeutrinoAM, mnatnbE, true},
	{ENeutrinoAM, enbamnE, false},

	{ENeutrinoT, tnatnE, false},
	{ENeutrinoT, enbatnE, true},
	{ENeutrinoT, mnbatnE, true},
	{ENeutrinoAT, tnatnE, false},
	{ENeutrinoAT, enbatnE, false},
	{ENeutrinoAT, mnbatnE, false},

	{EEndParticle,end_tE,false}
};

void WeakCoupling::Channel_T::Coef(CMatrixAddOnlyView& aCoef) const
{
	const Medium& backgr = *Background();
	double neutrinoConc =  backgr.neutrinoConc*backgr.neutrinoClusteringModifier;
	const int nn = Ranges().nE();
	const CBinning& E = Ranges().midE();
	double constMult = (BC().ss2-BC().ss_2)*neutrinoConc;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		double mult = constMult / E[iPrim];
		for(int iSec = 0; iSec<=iPrim; iSec++)
		{
			aCoef.Add(iSec, iPrim, mult*E[iSec]*fChannel.Sigma(iPrim, nn-1-(iPrim-iSec), fEffSecondary));
		}
	}
}

}//namespace couplings{
