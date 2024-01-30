#pragma once
#include "Vector.h"
#include "TParticle.h"
#include "ParticleList.h"
#include "Parameters.h"

class Medium;
class CouplingChannel;

enum RateMode
{
	RMdefault=0,
	RMtau,
	RMgrammage,

	RMEnd
};

class CouplingParameters : protected Parameters
{
public:
	inline bool Interacts(TParticle aParticle) const { return (bool)fInteractingParticles[aParticle]; }
	inline bool CoeffTestEnabled() { return COEF_TEST_ON; }
	CouplingParameters();
protected:
	static bool ready;

	static bool fInteractingParticles[EEndParticle];
	static bool fEnableCustomCouplings;
	static bool NO_PP; //do not take into account Par Production by photons (PP)
	static bool NO_ICS; //do not take into account Inverse Compton Scatering by electrons (ICS)
	static bool NO_TPP; //do not take into account Triple Par Production by electrons (TPP)
	static bool NO_DPP; //do not take into account Double Pair Production by photons (DPP)
	static bool NO_SYNCHROTRON; //do not take into account synchrotron radiation of electrons (effect only on photons)
	static bool NO_PHOTON_SPLITTING;//disable photon splitting (process violating Lorentz invariance)
	static bool NO_PPP; // do not take into account Par Production by Protons (PPP)
	static bool NO_PI; //do not take into account pion production
	static bool NO_N_DECAY;  // do not take into account neutron decay
	static bool NO_PHOTODISINTEGRATION;//do not take into account photodisintegration of nuclei
	static bool NO_pA;//do not take into account interactions of nuclei with proton gas (experimental)
	static bool NO_pp;//do not take into account interactions of protons with proton gas
	static bool JOIN_ELECTRON_POSITRON;//mix positrons as electrons (decreases number of dimensions in Jacobian, positrons should be switched off)
	static bool NO_MU_PP;//turn off muon pair production by photons
	static RateMode LOCAL_pp;
	static std::string LOCAL_pp_rate;//table with local pp tau (for use with 1Mpc short distance test)
	static std::string pA_grammage;

	static bool NO_NEUTRINO_t;//turn off neutrino interactions t-chanel
	static bool NO_NEUTRINO_s_HADRON;//turn off neutrino interactions s-chanel (hadron part)
	static bool NO_NEUTRINO_s_NON_HADRON;//turn of neutrino interactions s-chanel (nonhadron part)
	static bool NO_NEUTRINO_s;//turn off neutrino interactions s-chanel (non-hadron part)
	static bool COEF_TEST_ON;
};



class Coupling
{
public:
	Coupling(void);
	//to do: move all coefficients from PropagCoef to Coupling
	// PropagCoef will be used to store background information only (magnetic field, photon and neutrino backgrounds)
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	//SetBackgrounds should be called prior to calling DebugOutput()
	virtual void DebugOutput(const char* aFolder) const {};
	virtual ~Coupling(void);
	inline const CPointerArray<CouplingChannel>&	Channels() const {return fChannels;}
	bool IsEnabled() const;
	inline const Medium* Background() const {return fBackground;}
protected:
	inline void AddChannel(CouplingChannel* aChannel) {fChannels.add(aChannel);}
	const Medium*						fBackground;
private:
	CAutoDeletePtrArray<CouplingChannel>		fChannels;
};

class CouplingChannel
{
protected:
	CouplingChannel(TParticle aPrimary, TParticle aSecondary):
		 fPrimary(aPrimary),fSecondary(aSecondary){};
public:
	virtual ~CouplingChannel(){}
	inline TParticle Primary() const {return fPrimary; }
	inline TParticle Secondary() const {return fSecondary; }
	inline bool IsEnabled() const;
	virtual void Coef(CMatrixAddOnlyView& aCoef) const = 0;
protected:
	TParticle fPrimary;
	TParticle fSecondary;
};

template<class TCoupling>
class CouplingChannelT : public CouplingChannel
{
protected:
	CouplingChannelT(TCoupling* aCoupling, TParticle aPrimary, TParticle aSecondary):
	CouplingChannel(aPrimary, aSecondary),fCoupling(*aCoupling){};
public:
	virtual void Coef(CMatrixAddOnlyView& aCoef) const = 0;
protected:
	TCoupling& fCoupling;
	inline const Medium* Background() const { return fCoupling.Background();}
};

class CouplingList : public CouplingParameters
{
public:
	CouplingList();
	virtual ~CouplingList();

	inline static CouplingList* Instance(){
		if(!fInstance)
			fInstance = new CouplingList();
		return fInstance;
	};

	inline static CouplingList* Set(CouplingList* aCouplingList){
		CouplingList* result = fInstance;
		fInstance = aCouplingList;
		return result;
	};

	inline CPointerArray<Coupling>&	Couplings() {return fCouplings;}
	void DebugOutput(const char* aFolder) const;
	void init();

protected:
	//defined in CustomCouplings.cpp, called from init()
	void CustomInit();
	static CouplingList*					fInstance;
	CAutoDeletePtrArray<Coupling>			fCouplings;
};

inline bool CouplingChannel::IsEnabled() const
{
	CParticleList* pf = CParticleList::Instance();
	CouplingList& cl = *CouplingList::Instance();
	return cl.Interacts(fPrimary) && pf->IsEnabled(fPrimary) && pf->IsEnabled(fSecondary);
}
