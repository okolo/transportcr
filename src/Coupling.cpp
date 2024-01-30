#include "Coupling.h"
#include "Parameters.h"
#include "ICS.h"
#include "PP.h"
#include "Synchrotron.h"
#include "TPP.h"
#include "Dpp.h"
#include "MPP.h"
#include "Ppp.h"
#include "GZK2.h"
#include "PhotoDisintegration.h"
#include "NeutronDecay.h"
#include "Weak.h"
#include "AA.h"
#include "AAold.h"
//#include "Deflection.h"
#include "GammaSplitting.h"

using namespace couplings;

Coupling::Coupling(void)
{
}

Coupling::~Coupling(void)
{
}

bool Coupling::IsEnabled() const
{
	for(int i=fChannels.length()-1; i>=0; i--)
		if(fChannels(i).IsEnabled())
			return true;
	return false;
}

void Coupling::SetBackgrounds(const Medium& aPropagCoef)
{
	fBackground = &aPropagCoef;
}

CouplingList::CouplingList()
{

}

CouplingList::~CouplingList()
{

}

enum TModeICS
{
	ICSModeOld = 0,
	ICSModeCEL,
	ICSModeCELsplitted,

	ICSModeEnd //last element
};

enum TModePP
{
	PPModeOld = 0,
	PPModeApril2010,

	PPModeEnd //last element
};


CouplingList* CouplingList::fInstance;

void CouplingList::init()
{
	TModeICS ICS_Mode = ICSModeCEL; //other modes left only for comparison with old results
	TModePP PP_Mode = PPModeApril2010; //other modes left only for comparison with old results

	if(!NO_PP)
		switch(PP_Mode)
		{
		case PPModeOld:
			fCouplings.add(new PPold());
			break;
		case PPModeApril2010:
			fCouplings.add(new PPapr10());
			break;
		default:
			ThrowError("Unexpected PP_Mode value " + ToString(ICS_Mode));
		}
	if(!NO_ICS)
		switch(ICS_Mode)
		{
		case ICSModeOld:
			fCouplings.add(new ICSold());
			break;
		case ICSModeCEL:
			fCouplings.add(new ICScel());
			break;
		case ICSModeCELsplitted:
			fCouplings.add(new ICSsplitted());
			break;
		default:
			ThrowError("Unexpected ICS_Mode value " + ToString(ICS_Mode));
		}
	if(!NO_TPP)
		fCouplings.add(new TPP());
	if(!NO_DPP)
		fCouplings.add(new Dpp());
	if(!NO_SYNCHROTRON)
		fCouplings.add(new Synchrotron());
	if(!NO_MU_PP)
		fCouplings.add(new MPP());
	if(!NO_PPP)
		fCouplings.add(new PPP());
	if(!NO_PI)
		fCouplings.add(new GZK());
	if(!NO_N_DECAY)
		fCouplings.add(new NeutronDecay());
	if(!NO_PHOTODISINTEGRATION)
		fCouplings.add(new PhotoDisintegration());
	if(!(NO_NEUTRINO_s && NO_NEUTRINO_t))
		fCouplings.add(new WeakCoupling());
    if(!NO_pA)
		fCouplings.add(new NucleiProton(pA_grammage.c_str()));
    if(!NO_pp)
    	fCouplings.add(new NucleonProton());
    if(LOCAL_pp != RMdefault)
    	fCouplings.add(new NucleonProton(LOCAL_pp, PPCodeKachelriess, LOCAL_pp_rate.c_str()));
    if(!NO_PHOTON_SPLITTING)
        fCouplings.add(new GammaSplitting());
    if(fEnableCustomCouplings)
    	CustomInit();
}

void CouplingList::DebugOutput(const char* aFolder) const
{
	for(int i=0; i<fCouplings.length(); i++)
		if(fCouplings(i).IsEnabled())
			fCouplings(i).DebugOutput(aFolder);
}

bool CouplingParameters::ready = false;
bool CouplingParameters::fInteractingParticles[EEndParticle];
bool CouplingParameters::fEnableCustomCouplings = false;
bool CouplingParameters::NO_PP = false; //do not take into account Par Production by photons (PP)
bool CouplingParameters::NO_ICS = false; //do not take into account Inverse Compton Scatering by electrons (ICS)
bool CouplingParameters::NO_TPP = false; //do not take into account Triple Par Production by electrons (TPP)
bool CouplingParameters::NO_DPP = false; //do not take into account Double Pair Production by photons (DPP)
bool CouplingParameters::NO_SYNCHROTRON = false; //do not take into account synchrotron radiation of electrons (effect only on photons)
bool CouplingParameters::NO_PHOTON_SPLITTING = true;//disable photon splitting (process violating Lorentz invariance)
bool CouplingParameters::NO_PPP = false; // do not take into account Par Production by Protons (PPP)
bool CouplingParameters::NO_PI = false; //do not take into account pion production
bool CouplingParameters::NO_N_DECAY = false;  // do not take into account neutron decay

bool CouplingParameters::NO_PHOTODISINTEGRATION = false;//do not take into account photodisintegration of nuclei
bool CouplingParameters::NO_pA = true;//do not take into account interactions of nuclei with proton gas (experimental)
bool CouplingParameters::NO_pp = true;//do not take into account interactions of protons with proton gas
bool CouplingParameters::JOIN_ELECTRON_POSITRON = false;//mix positrons as electrons (decreases number of dimensions in Jacobian, positrons should be switched off)
bool CouplingParameters::NO_MU_PP = true;//turn off muon pair production by photons

bool CouplingParameters::NO_NEUTRINO_t = true;//turn off neutrino interactions t-chanel
bool CouplingParameters::NO_NEUTRINO_s_HADRON = true;//turn off neutrino interactions s-chanel (hadron part)
bool CouplingParameters::NO_NEUTRINO_s_NON_HADRON = true;//turn of neutrino interactions s-chanel (nonhadron part)
bool CouplingParameters::NO_NEUTRINO_s = true;//turn off neutrino interactions s-chanel (non-hadron part)
bool CouplingParameters::COEF_TEST_ON = false;//output debug coupling information and background
RateMode CouplingParameters::LOCAL_pp = RMdefault;//local pp interaction with rate given by LOCAL_pp_tau table in Mpc^-1
std::string CouplingParameters::LOCAL_pp_rate = "local_pp";//table with local pp rate (in Mpc^-1 or if tau is given 1Mpc short distance test should be used)
std::string CouplingParameters::pA_grammage = "";

CouplingParameters::CouplingParameters()
{
	if(ready) return;

	for(int i=0; i<EEndParticle; i++)
		fInteractingParticles[i] = true;

	READ_BOOL_SETTING(NO_PP);
	READ_BOOL_SETTING(NO_ICS);
	READ_BOOL_SETTING(NO_TPP);
	READ_BOOL_SETTING(NO_DPP);
	READ_BOOL_SETTING(NO_SYNCHROTRON);
    READ_BOOL_SETTING(NO_PHOTON_SPLITTING);

	READ_BOOL_SETTING(NO_PPP);
	READ_BOOL_SETTING(NO_MU_PP);
	READ_BOOL_SETTING(NO_PI);
	READ_BOOL_SETTING(NO_N_DECAY);

	READ_BOOL_SETTING(NO_PHOTODISINTEGRATION);
	READ_BOOL_SETTING(NO_pp);
    READ_BOOL_SETTING(NO_pA);
    if(!NO_pA) {
		if(LOCAL_pp || (!NO_pp))
			ThrowError("Either pp or pA interactions can be enabled, but not both (pA icludes pp");
		READ_STRING_SETTING(pA_grammage);
		if(pA_grammage.length()>0){
			CopyFileOrDir(pA_grammage, plt_local_dir);
		}
	}
    READ_SWITCH_SETTING(LOCAL_pp, RMEnd);
    if(LOCAL_pp!=RMdefault)
    {
    	READ_STRING_SETTING(LOCAL_pp_rate);
    	CopyFileOrDir(LOCAL_pp_rate, plt_local_dir);
    }

	READ_BOOL_SETTING(COEF_TEST_ON);

	READ_BOOL_SETTING(NO_NEUTRINO_t);
	READ_BOOL_SETTING(NO_NEUTRINO_s_HADRON);
	READ_BOOL_SETTING(NO_NEUTRINO_s_NON_HADRON);
	NO_NEUTRINO_s = NO_NEUTRINO_s_HADRON&&NO_NEUTRINO_s_NON_HADRON;
	
	READ_BOOL_SETTING(JOIN_ELECTRON_POSITRON);

	Reader()->readBoolPar("EnableCustomCouplings", fEnableCustomCouplings);
	for(int i=0; i<EEndParticle; i++)
	{
		std::string key = ToString("int_") + ParticleData::getParticleFileName((TParticle)i);
		Reader()->readBoolPar( key.c_str(), fInteractingParticles[i]);
	}
	ready = true;
}
