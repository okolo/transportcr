// MassiveNeutrino.h: 
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MASSIVENEUTRINO_H__6DF110E3_4AEB_11D5_B6B8_E8293AD9338A__INCLUDED_)
#define AFX_MASSIVENEUTRINO_H__6DF110E3_4AEB_11D5_B6B8_E8293AD9338A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "LEPData.h"
#include "Sigma.h"
#include "SecondaryDifSigma.h"
#include "Ranges.h"
#include "Parameters.h"

class NuSigmaParameters : protected Parameters
{
protected:
	NuSigmaParameters();
	static double NEUTRINO_S_CHANNEL_MULTIPLIER;//s-channel cross section multiplier, need for debugging purposes
	static double NEUTRINO_T_CHANNEL_MULTIPLIER;//t-channel cross section multiplier, need for debugging purposes
	static int RESONANCE_RESOLUTION_ACCURACY;
	static bool Z_CHANEL_YOSHIDA;//calculate neutrino interactions s-chanel using Yoshida's (wrong) formula (astro-ph/9608186)
	static double m_en;//electron neutrino mass in internal units
	static double m_mn;//muon neutrino mass in internal units
	static double m_tn;//tau neutrino mass in internal units
};

class CNeutrinoSChanelSigma : public CSigma, public CDifSigma, NuSigmaParameters
{
public:
	double SigmaTotalResonanceValue();
	CNeutrinoSChanelSigma(double gL, double gR);
	//from ISigma
	double Sigma(double _s) const;
	//CDifSigma
	double Sigma(double _s, double _ratio) const;

	/** Calculate (Integral(sigma(S)dS)/deltaS), that is mean sigma in the resonance
	energy bin, the integral is then stored in CNeutrinoSChanelSigma::resonanceIntegral variable.
	Returns resonance energy bin number
	*/
	int CalculateResonanceIntegrals(double sMin);
private:

// mean sigma in resonance bin	(Integral(sigma(S)dS)/deltaS)
	double resonanceIntegral;

// mean sigma in neighbour bins
	double leftResonanceIntegral;
	double rightResonanceIntegral;
//resonance bin bounds
	double sMinResonance;
	double sMaxResonance;
//neighbour bins bounds
	double sLeftResonance;//=sMinResonance/BC().ss1
	double sRightResonance;//=sMaxResonance*BC().ss1

	double m_gL2;//gL^2
	double m_gR2;//gR^2
	double m_C1;//==12*Pi/Mz^2*(GammaIn/GammaZ)*GammaZ^2
	double m_C2;//==3/2/(gL^2+gR^2)

protected:
	static const double s_mult; // multiplier obtained from comparison of right
	// formula (see Gondolo & Co. hep-ph/9209236) with one taken from Yoshida's
	// paper astro-ph/9608186
};

class CZHadronChanel
{
public:
	CZHadronChanel();
	static void InitLEPData();
	static void CloseLEPData();
protected:
	static CLEPData* m_LEPData;
};

class CNeutrinoTChanelSigma : public CSecondaryDifSigma, public CSigma, public CDifSigma, NuSigmaParameters
{
public:
	CNeutrinoTChanelSigma(double _Smin, int _nn=Ranges().nE());
	virtual ~CNeutrinoTChanelSigma();
	//from ISigma
	double Sigma(double _s) const;
	//CDifSigma
	double Sigma(double _s,double _ratio) const;
	
//sum of Sigmas on all background neutrinos at given UHECR energy (massive neutrino case)
	inline static double SigmaTotMassive(int _E){ASSERT(_E<Ranges().nE()&&m_SigmaTotMassive);return m_SigmaTotMassive[_E];}
	double SigmaTotMassive(double _E) const;
//from CSecondaryDifSigma
	virtual void Init();
protected:
	double SigmaLeading(double _s,double _ratio) const;
	static MuDecaySpectrum m_muDecay;
	static double* m_SigmaTotMassive;
//	static MuDecaySpectrum m_tauDecay;//temporay, later MuDecaySpectrum should be
	//replaced by TauDecaySpectrum
};

class CMassiveNeutrinoSHadronChanel : public CNeutrinoSChanelSigma,
public CZHadronChanel
{
public:
	CMassiveNeutrinoSHadronChanel();
	virtual ~CMassiveNeutrinoSHadronChanel();
	void Init(double _mBackgr);
	double P(double _E, double _nE, double _mBackgr, CLEPData::EParticles _secParticle) const;
	inline double R(double _E, double _mBackgr) const;
	inline double P(int _E, int _nE, CLEPData::EParticles _secParticle) const;
	inline double R(int _E) const;

protected:
	double PMultiplier(double _E, double _nE, double _mBackgr, CLEPData::EParticles _secParticle) const;
	bool isReady;
	matrix* m_P;
	PDouble m_R;
};

inline double CMassiveNeutrinoSHadronChanel::R(double _E, double _mBackgr) const
{
	return Sigma(2.*_E*_mBackgr);
}

inline double CMassiveNeutrinoSHadronChanel::P(int _E, int _nE, CLEPData::EParticles _secParticle) const
{
	ASSERT(_nE<=_E);
	ASSERT(_E<Ranges().nE());
	ASSERT(isReady);
	return m_P[_secParticle][_E][_nE];
}

inline double CMassiveNeutrinoSHadronChanel::R(int _E) const
{
	ASSERT(_E<Ranges().nE());
	ASSERT(isReady);
	return m_R[_E];
}


class CNeutrinoSNonHadronChanel
{
public:
	CNeutrinoSNonHadronChanel(double _Smin, int _nn=Ranges().nE());
	virtual void Init();//set link from muon secondary chanel to electron one
	inline double Sigma(int _s,int _ratio,TParticle _particle) const;//sigma multiplied by E (initial particle energy)
	virtual ~CNeutrinoSNonHadronChanel(){};
protected:
//overrided from CSecondaryDifSigma
	CNeutrinoSChanelSigma m_neutrinoSChanelSigma;
	CNeutrinoSChanelSigma m_leptonSChanelSigma;
	MuDecaySpectrum m_muDecay;
	CSecondaryDifSigma m_leptonSChanel;
	CSecondaryDifSigma m_neutrinoSChanel;
};

inline double CNeutrinoSNonHadronChanel::Sigma(int _s,int _ratio,TParticle _particle) const
{
	return (m_leptonSChanel.Sigma(_s, _ratio,_particle)+m_neutrinoSChanel.Sigma(_s, _ratio,_particle));
};

class CEnAenChanel : public CNeutrinoTChanelSigma
{
public:
	CEnAenChanel(double _Smin, int _nn=Ranges().nE());
};

#ifndef NEUTRINO_T_PRIMARY_CHANNELS_ONLY

class CEnAmnChanel : public CNeutrinoTChanelSigma
{
public:
	CEnAmnChanel(double _Smin, int _nn=Ranges().nE());
protected:
	void PostInit();
};

class CMnAmnChanel : public CNeutrinoTChanelSigma
{
public:
	CMnAmnChanel(double _Smin, int _nn=Ranges().nE());
protected:
	void PostInit();
};

#else

typedef CEnAenChanel CEnAmnChanel;
typedef CEnAenChanel CMnAmnChanel;


#endif	//NEUTRINO_T_PRIMARY_CHANNELS_ONLY

//since currently we use the same decay spectrum for mu and tau
//we can use the same classes for the chanels with tau neutrino:
typedef CMnAmnChanel CMnAtnChanel;
typedef CEnAmnChanel CEnAtnChanel;
typedef CMnAmnChanel CTnAtnChanel;

#endif // !defined(AFX_MASSIVENEUTRINO_H__6DF110E3_4AEB_11D5_B6B8_E8293AD9338A__INCLUDED_)
