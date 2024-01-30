#pragma once
#include "Coupling.h"
#include "LEPData.h"
#include "MassiveNeutrino.h"

extern int NO_NEUTRINO_s_HADRON;
extern int NO_NEUTRINO_s_NON_HADRON;

class CSecondaryDifSigma;        
class CNeutrinoTChanelSigma;         
class CNeutrinoSNonHadronChanel;    
class CMassiveNeutrinoSHadronChanel;

namespace couplings{

class WeakCoupling : public Coupling, NuSigmaParameters, CouplingParameters
{
	enum TNeutrinoTChanels
	{
		enaenE=0,
		enamnbE,
		enbamnE,
		enatnbE,
		mnamnE,
		mnatnbE,

		end_tE,

		mnbatnE = mnamnE,
		tnatnE = mnatnbE,
		enbatnE = enbamnE
	};

	struct TChannelMapEntry
	{
		TParticle			Primary;
		TNeutrinoTChanels	Channel;
		bool				UseAntiparticleCoef;
	};

	enum TNeutrinoSChanels
	{
		s_eE=0,
		s_mE,
		s_tE,

		end_sE
	};

	class Channel_AbsorptionS : public CouplingChannelT<WeakCoupling>
	{
	public:
		Channel_AbsorptionS(WeakCoupling* aCoupling, TParticle aPrim);
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		const CMassiveNeutrinoSHadronChanel&	fChannel;
		double									fMult;	
	};

	class Channel_AbsorptionT : public CouplingChannelT<WeakCoupling>
	{
	public:
		Channel_AbsorptionT(WeakCoupling* aCoupling, TParticle aPrim):
		  CouplingChannelT<WeakCoupling>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;	
	};

	class Channel_HadronS : public CouplingChannelT<WeakCoupling>
	{
	public:
		Channel_HadronS(WeakCoupling* aCoupling, TParticle aPrim, TParticle aSec);
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		const CMassiveNeutrinoSHadronChanel& fChannel;
		CLEPData::EParticles	fSecondaryGroup;
		double					fMult;
	};

	class Channel_NonHadronS : public CouplingChannelT<WeakCoupling>
	{
	public:
		Channel_NonHadronS(WeakCoupling* aCoupling, TParticle aPrim, TParticle aSec);
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		const CNeutrinoSNonHadronChanel& fChannel;
	};

	class Channel_T : public CouplingChannelT<WeakCoupling>
	{
	public:
		Channel_T(WeakCoupling* aCoupling, const TChannelMapEntry& aEntry, TParticle aSec);
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		const CSecondaryDifSigma&	fChannel;
		TParticle					fEffSecondary;
	};

public:
	WeakCoupling(void);
	virtual ~WeakCoupling();
private:
	void Init();
	void Close();
	void WriteCoef();
	static TNeutrinoSChanels GetSChanelIndex(TParticle aPrimary);

	CSecondaryDifSigma**           m_neutrinoChanals;
	CNeutrinoTChanelSigma*         m_tChanelSigma;
	CNeutrinoSNonHadronChanel**    m_sNonHadronChanels;
	CMassiveNeutrinoSHadronChanel* m_MassiveNeutrinoSHadronChanels;
	static const TChannelMapEntry fTChannelMap[];
};

}//namespace couplings
