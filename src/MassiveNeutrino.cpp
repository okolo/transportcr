#include "MassiveNeutrino.h"
#include "const.h"
#include <math.h>
#include "Units.h"


double NuSigmaParameters::NEUTRINO_S_CHANNEL_MULTIPLIER = 1.;//s-channel cross section multiplier, need for debugging purposes
double NuSigmaParameters::NEUTRINO_T_CHANNEL_MULTIPLIER = 1.;//t-channel cross section multiplier, need for debugging purposes
bool NuSigmaParameters::Z_CHANEL_YOSHIDA = 0;//calculate neutrino interactions s-chanel using Yoshida's (wrong) formula (astro-ph/9608186)
int NuSigmaParameters::RESONANCE_RESOLUTION_ACCURACY = 1000;
double NuSigmaParameters::m_en=-1;//electron neutrino mass in internal units
double NuSigmaParameters::m_mn=-1;//muon neutrino mass in internal units
double NuSigmaParameters::m_tn=-1;//tau neutrino mass in internal units

NuSigmaParameters::NuSigmaParameters()
{
	if(m_en>0)
		return;//already initialized

	//default values in eV
	m_en = 1;
	m_mn = 1;
	m_tn = 1;

	READ_BOOL_SETTING(Z_CHANEL_YOSHIDA);
	READ_DOUBLE_SETTING(NEUTRINO_S_CHANNEL_MULTIPLIER);
	READ_DOUBLE_SETTING(NEUTRINO_T_CHANNEL_MULTIPLIER);
	READ_INT_SETTING(RESONANCE_RESOLUTION_ACCURACY);

	READ_DOUBLE_SETTING(m_en);
	READ_DOUBLE_SETTING(m_mn);
	READ_DOUBLE_SETTING(m_tn);

	//converting to internal units
	m_en *= (1e-6/units.Eunit);
	m_mn *= (1e-6/units.Eunit);
	m_tn *= (1e-6/units.Eunit);
}

CNeutrinoSChanelSigma::CNeutrinoSChanelSigma(double gL, double gR)
:resonanceIntegral(-1.),
sMinResonance(0),
sMaxResonance(0)
{
extern const double Pi;
	m_gL2 = gL*gL;
	m_gR2 = gR*gR;
	m_C1 = 12.*Pi/FWeak::Mz2*FWeak::gIn_en*FWeak::GammaZ2*s_mult;
	m_C2 = 3./(m_gL2+m_gR2);//3/(gL^2+gR^2)
}
const double CNeutrinoSChanelSigma::s_mult = 3.87;

double CNeutrinoSChanelSigma::Sigma(double _s) const
//multiplied by Gamma_in * Gamma_z (that is without Gamma_out/Gamma_z)
{
	if((resonanceIntegral>0)&&(_s>sLeftResonance)&&(_s<sRightResonance))
	{
		if(_s<sMinResonance)
			return leftResonanceIntegral;
		if(_s>sMaxResonance)
			return rightResonanceIntegral;
		return resonanceIntegral;
	}
	double delta=_s-FWeak::Mz2;
	double result=m_C1*NEUTRINO_S_CHANNEL_MULTIPLIER;
	if(Z_CHANEL_YOSHIDA)
		result*=(_s/(delta*delta+FWeak::Cz0));
	else
		result*=(_s/(delta*delta+_s*_s*FWeak::Cz5));
	ASSERT(result>=0.0);
	return result;
}

double CNeutrinoSChanelSigma::SigmaTotalResonanceValue()
{
// switching of resonance integral approximation
	double mem = resonanceIntegral;
	resonanceIntegral = -1.;

	double result = Sigma(FWeak::Mz2);
	
// switching on resonance integral approximation	
	resonanceIntegral = mem;
	
	return result;
}


double CNeutrinoSChanelSigma::Sigma(double _s,double _ratio) const
//multiplied by Gamma_in * Gamma_z (that is without Gamma_out/Gamma_z)
{
	ASSERT(_ratio<=1.);
	double r1=1.-_ratio;
	double sum = r1*r1 + _ratio*_ratio;
	//double result = m_C2*(m_gL2*_ratio*_ratio + m_gR2*r1*r1);//this is dif cross section for leading particle
	double result = 0.5*m_C2*sum*(m_gL2 + m_gR2);//this averaged dif cross section == 1/2 * (Sigma_leading + Sigma_nonleading)
	result*=Sigma(_s);
	ASSERT(result>=0.0);
	return result;
}

int RESONANCE_RESOLUTION_ACCURACY = 1000;

int CNeutrinoSChanelSigma::CalculateResonanceIntegrals(double sMin)
{
	resonanceIntegral = -1.;//invalidating previous resonance integration if any
	int binNo;
	{
		const double s = BC().s;
//calculating resonance bin number
		binNo = (int)(s*log10(FWeak::Mz2/sMin));

//calculating integral limits
		sMinResonance = sMin*pow(10.,binNo/s);
		sMaxResonance = sMinResonance*BC().ss1;
		sLeftResonance = sMinResonance/BC().ss1;
		sRightResonance = sMaxResonance*BC().ss1;
	}
	if(RESONANCE_RESOLUTION_ACCURACY<=0)//resonance integration was turned off explicitly
	{
		return binNo;
	}
	leftResonanceIntegral = MeanSigma(sLeftResonance,sMinResonance,RESONANCE_RESOLUTION_ACCURACY);
	rightResonanceIntegral = MeanSigma(sMaxResonance,sRightResonance,RESONANCE_RESOLUTION_ACCURACY);
	
	//must be last (important)
	resonanceIntegral = MeanSigma(sMinResonance,sMaxResonance,RESONANCE_RESOLUTION_ACCURACY);
		
	return binNo;
}


double CNeutrinoTChanelSigma::Sigma(double _s) const
//total cross section for t-chanel
{
	double result=0;
#define M2 FWeak::Mw2
	double x = _s/M2;
	if(x>1e-3)
		result = 2.*(1.+0.5*x-(1./x+1.)*log(1.+x));
	else
		result = 0.33333333*x*x;
#undef M2
	ASSERT(result>=0.);
	return result*FWeak::Cw2/_s*NEUTRINO_T_CHANNEL_MULTIPLIER;
}

double* CNeutrinoTChanelSigma::m_SigmaTotMassive = NULL;

CNeutrinoTChanelSigma::CNeutrinoTChanelSigma(double _Smin, int _nn)
:CSecondaryDifSigma(NULL,_Smin,_nn)
{
	m_difSigma = this;
}

CNeutrinoTChanelSigma::~CNeutrinoTChanelSigma()
{
	delete[] m_SigmaTotMassive;
	m_SigmaTotMassive = NULL;
}

double CNeutrinoTChanelSigma::Sigma(double _s,double _ratio) const
{
	if(!this)
		return 0.;//interactions were turned off
	double result=0.5*(SigmaLeading(_s,_ratio)+SigmaLeading(_s,1.-_ratio));
	//double result=SigmaLeading(_s,_ratio);
	result *=(FWeak::Cw2*_s*NEUTRINO_T_CHANNEL_MULTIPLIER);
	ASSERT_VALID_NO(result);
	return result;
}

//Leading particle sigma divided by (FWeak::Cw2*s)
//SigmaNonLeading(s,ratio) = SigmaLeading(s,1-ratio)
double CNeutrinoTChanelSigma::SigmaLeading(double _s,double _ratio) const
{
	double divisor=_s*(1.-_ratio)+FWeak::Mw2;
	divisor*=divisor;
	double result=_ratio*_ratio/divisor;
	ASSERT_VALID_NO(result);
	return result;
}

double CNeutrinoTChanelSigma::SigmaTotMassive(double _E) const
{
	double result = Sigma(2.*m_en*_E);
	result += Sigma(2.*m_mn*_E);
	result += Sigma(2.*m_tn*_E);
	return result;
}

void CNeutrinoTChanelSigma::Init()
{
	const double Emin = Ranges().Emin();
	const int nn = Ranges().nE();
	CSecondaryDifSigma::Init();
	
	if(m_SigmaTotMassive)
		return;
	
	m_SigmaTotMassive = new double[nn];
	double E=Emin*BC().ss2;
	int i;
	for(i=0;i<nn;i++,E*=BC().ss1)
	{
		m_SigmaTotMassive[i] = SigmaTotMassive(E);
	}
/*	printSigmaTot(PLT_DIR"dif", m_Smin, nn, 100000);
	Print(PLT_DIR"tot", m_Smin, nn);*/
}

CMassiveNeutrinoSHadronChanel::CMassiveNeutrinoSHadronChanel():
CNeutrinoSChanelSigma(1.,1.),//in hadron chanel gL & gR are not used
isReady(false),
m_P(NULL),
m_R(NULL)
{
}

CMassiveNeutrinoSHadronChanel::~CMassiveNeutrinoSHadronChanel()
{
	const int nn = Ranges().nE();
	if(m_P)
		for(int i=0;i<CLEPData::endE;i++)
		{
			for(int j=0;j<nn;j++)
				delete[] m_P[i][j];
			delete m_P[i];
		}
		
	delete m_P;
	delete[] m_R;
}

double CMassiveNeutrinoSHadronChanel::P(double _E, double _nE, double _mBackgr, CLEPData::EParticles _secParticle) const
{
	double result = R(_E,_mBackgr)*PMultiplier(_E, _nE, _mBackgr, _secParticle);//Sigma replaced by R (factor 2)
	return result;
}


double CMassiveNeutrinoSHadronChanel::PMultiplier(double _E, double _nE, double _mBackgr, CLEPData::EParticles _secParticle) const
{
	double gamma = sqrt(0.5*_E/_mBackgr);
	double result = FWeak::gOut_q*m_LEPData->N(_secParticle,_nE/FWeak::Mz/gamma)/FWeak::Mz/gamma;
	return result;
}


void CMassiveNeutrinoSHadronChanel::Init(double _mBackgr)
{
	const double Emin = Ranges().Emin();
	const int nn = Ranges().nE();
	if(isReady)
		return;
	//int resIndex = CalculateResonanceIntegrals(2.*Emin*_mBackgr);
	m_P = new matrix[CLEPData::endE];
	m_R = new double[nn];
	int i;
	for(i=0;i<CLEPData::endE;i++)
	{
		m_P[i]=new PDouble[nn];
		for(int j=0;j<nn;j++)
			m_P[i][j]=new double[j+1];
	}
	double E = Emin*BC().ss2;
	for(int j=0;j<nn;j++,E*=BC().ss1)
	{
		m_R[j]=R(E,_mBackgr);
		double nE = Emin*BC().ss2;
		for(int k=0;k<=j;k++,nE*=BC().ss1)
			for(i=0;i<CLEPData::endE;i++)
				m_P[i][j][k] = m_R[j] * PMultiplier(E, nE, _mBackgr, (CLEPData::EParticles)i);
	}
	isReady = true;
}

void CZHadronChanel::InitLEPData()
{
	if(m_LEPData!=NULL) return;
	m_LEPData = new CLEPData;
	m_LEPData->read();
}

void CZHadronChanel::CloseLEPData()
{
	delete m_LEPData;
	m_LEPData=NULL;
}

CZHadronChanel::CZHadronChanel()
{
	InitLEPData();
}

CLEPData* CZHadronChanel::m_LEPData = NULL;

CNeutrinoSNonHadronChanel::CNeutrinoSNonHadronChanel(double _Smin, int _nn):
m_neutrinoSChanelSigma(FWeak::gL_n,FWeak::gR_n),
m_leptonSChanelSigma(FWeak::gL_l,FWeak::gR_l),
m_leptonSChanel(&m_leptonSChanelSigma,_Smin,_nn),
m_neutrinoSChanel(&m_neutrinoSChanelSigma,_Smin,_nn)
{
	m_neutrinoSChanelSigma.CalculateResonanceIntegrals(_Smin);
	m_leptonSChanelSigma.CalculateResonanceIntegrals(_Smin);

	//m_difSigma = &m_neutrinoSChanelSigma;
	m_neutrinoSChanel.AddPrimaryMode(ENeutrinoE,FWeak::gIn_en);
	m_neutrinoSChanel.AddPrimaryMode(ENeutrinoM,FWeak::gIn_mn);
	m_neutrinoSChanel.AddPrimaryMode(ENeutrinoT,FWeak::gIn_tn);


	m_leptonSChanel.AddPrimaryMode(EElectron,FWeak::gOut_e);

	m_leptonSChanel.AddDecayMode(&m_muDecay,FWeak::gOut_m+FWeak::gOut_t);//this is true only if CTauDecaySpectrum = CMuDecaySpectrum

	m_leptonSChanel.AddSecondaryParticle(EElectron);
	m_leptonSChanel.AddSecondaryParticle(ENeutrinoAE);
}

void CNeutrinoSNonHadronChanel::Init()//set link from muon secondary chanel to electron one
{
	m_leptonSChanel.Init();
	m_neutrinoSChanel.Init();

//	m_neutrinoSChanel.AddSecondaryLink(ENeutrinoM,EElectron);//improved on 21-st of Dec 2001
	m_leptonSChanel.AddSecondaryLink(ENeutrinoM,EElectron);//this is true only if CTauDecaySpectrum = CMuDecaySpectrum

	m_leptonSChanel.MoveSecondaryLink(ENeutrinoAE,ENeutrinoE);
	
	m_neutrinoSChanel.MakeAntiparticlesLikeParticles();
	m_leptonSChanel.MakeAntiparticlesLikeParticles();
//test
#if defined(_DEBUG)
/*	int nS = m_nn/2;
	double S = m_sMin*pow(BC().ss1,nS);
	double result = m_leptonChanel.CalculateSigma_tot(nS);
	double compareWith = 2*m_neutrinoSChanelSigma.Sigma(S);*/
	//((CSigma*)((ISigma*)(&m_neutrinoSChanelSigma)))->Print(PLT_DIR"sSigma",m_Smin,m_nn,2.);
#endif
}

/**
		t-chanels
*/

MuDecaySpectrum CNeutrinoTChanelSigma::m_muDecay;

#ifndef NEUTRINO_T_PRIMARY_CHANNELS_ONLY

CEnAmnChanel::CEnAmnChanel(double _Smin, int _nn)
:CNeutrinoTChanelSigma(_Smin, _nn)
{
	AddPrimaryMode(EElectron);
	AddDecayMode(&m_muDecay);
	AddSecondaryParticle(EElectron);
	AddSecondaryParticle(ENeutrinoAE);
}

void CEnAmnChanel::PostInit()
{
//	AddSecondaryLink(ENeutrinoAM,EElectron);//commented on 29.12.2001
	MoveSecondaryLink(EElectron,EPositron);
	MoveSecondaryLink(ENeutrinoAE,ENeutrinoE);
	AddSecondaryLink(ENeutrinoAM,EPositron);
#if defined(NU_DEBUG)
	int nS = m_nn/2;
	double S = m_Smin*pow(BC().ss1,nS);
	//double result = CalculateSigma_tot(nS);
	//double compareWith = Sigma(S);
	((CSigma*)((ISigma*)this))->Print((plt_local_dir + "tSigma").c_str(),m_Smin,m_nn);
#endif
}

CMnAmnChanel::CMnAmnChanel(double _Smin, int _nn)
:CNeutrinoTChanelSigma(_Smin, _nn)
{
	//AddDecayMode(&m_muDecay,2.);// 2 is not required since we use MakeAntiparticlesLikeParticles();
	AddDecayMode(&m_muDecay); //improved on 29.12.2001

	AddSecondaryParticle(EElectron);
	AddSecondaryParticle(ENeutrinoAE);
}

void CMnAmnChanel::PostInit()
{
	AddSecondaryLink(ENeutrinoM,EElectron);
	MoveSecondaryLink(ENeutrinoAE,ENeutrinoE);
	MakeAntiparticlesLikeParticles();
}

#endif //#ifndef NEUTRINO_T_PRIMARY_CHANNELS_ONLY

CEnAenChanel::CEnAenChanel(double _Smin, int _nn)
:CNeutrinoTChanelSigma(_Smin, _nn)
{
	AddPrimaryMode(EElectron);
	AddPrimaryMode(EPositron);
}


//#endif //NEUTRINO_DIMA
