// SecondaryDifSigma.cpp: implementation of the CSecondaryDifSigma class.
//
//////////////////////////////////////////////////////////////////////

#include "SecondaryDifSigma.h"
#include <math.h>
#include "ParticleList.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSecondaryDifSigma::CSecondaryDifSigma(CDifSigma* _primarySigma, double _Smin, int _nn)
:CDiscretDifSigma(_Smin,_nn),
m_difSigma(_primarySigma),
m_directRates(EEndParticle),
m_maxNoOfDecayModes(10),
m_noOfDecayModes(0),
isReady(0),
m_NoOfSigmasMem(0),
m_sigmas(EEndParticle),
m_sigmasMem(EEndParticle),
m_secondaryModes(NULL),
isPreLinkInitDone(0)
{
	const int nn = Ranges().nE();
	for(int i=0;i<EEndParticle;i++)
	{
		m_directRates[i]=0.;
//		m_sigmas[i]=0;//turn off particle (this turnes off both primary and all secondary chanels
//the appropriate particles must be turned on explisitly using AddSecondaryParticle() or/and AddPrimaryMode() function in derived classes
//		m_sigmasMem[i]=0;
	}
	m_secondaryModes = new TDecayMode[m_maxNoOfDecayModes];

	m_primarySigma.create(m_nn,nn);
}

CSecondaryDifSigma::~CSecondaryDifSigma()
{
	delete[] m_secondaryModes;
}

void CSecondaryDifSigma::AddDecayMode(DecaySpectrum *_mode, double _rate)
{
	ASSERT(m_maxNoOfDecayModes>m_noOfDecayModes);
	m_secondaryModes[m_noOfDecayModes].iRate = _rate;
	m_secondaryModes[m_noOfDecayModes].iSpectrum = _mode;
	m_noOfDecayModes++;
}

void CSecondaryDifSigma::AddPrimaryMode(TParticle _particle, double _rate)
{
	m_directRates[_particle] += _rate;
}


void CSecondaryDifSigma::AddSecondaryParticle(TParticle _particle)
{
	const int nn = Ranges().nE();
	if(m_sigmas[_particle]!=NULL)
		return;

	m_sigmasMem[m_NoOfSigmasMem] = new CMatrix(m_nn,nn);
/*	for(int i=0;i<m_nn;i++)
		m_sigmasMem[m_NoOfSigmasMem][i] = new double[nn];*/

	m_sigmas[_particle]=m_sigmasMem[m_NoOfSigmasMem];
	m_NoOfSigmasMem++;
}

void CSecondaryDifSigma::Init()
{
	int nS,nratio;
	double S,ratio,centerSigma,rightSigma,leftSigma;

	DECLARE_E_SHORTCUT

	const double ss1 = BC().ss1;
	const double ss2 = BC().ss2;
	const double ss_1 = BC().ss_1;
	const double s = BC().s;

	double integrationMultiplier = log(10.)/s/6.;


	for(nS = 0,S=m_Smin ; nS<m_nn; nS++,S*=ss1)
	{
		//calculating primary sigmas
		for(nratio = 0,ratio=Emin/Emax*ss1; nratio<nn; nratio++,ratio*=ss1)
			m_primarySigma[nS][nratio] = m_difSigma->Sigma(S,ratio);

		//calculating secondary sigmas
		CVector& primSigma = m_primarySigma[nS];
		for(nratio = 0,ratio=Emin/Emax*ss1; nratio<nn-1; nratio++,ratio*=ss1)
		{
			int nr1;
			double r1=ratio;
			double rLeft=1.;//=ratio/r1
			double rCenter=rLeft/ss2;//by def
			double rRight=rLeft/ss1;//by def
			for(nr1 = nratio; nr1<nn-1; nr1++,r1*=ss1,rLeft*=ss_1,rCenter*=ss_1,rRight*=ss_1)
			{
				leftSigma = primSigma[nr1];
				centerSigma = m_difSigma->Sigma(S,r1*ss2);
				rightSigma = primSigma[nr1+1];

				FOR_ALL_REAL_PARTICLES_INVOLVED(p)
				{
					//TParticle p = (TParticle)particle;
					if(!m_sigmas[p])
						continue;//skipping particles which were not added using AddSecondaryParticle
					double sum = 0;
					for(int decayMode=0;decayMode<m_noOfDecayModes;decayMode++)
					{
						double rate = m_secondaryModes[decayMode].iRate;
						DecaySpectrum* spec = m_secondaryModes[decayMode].iSpectrum;
						
						sum+=rate*( spec->N(rLeft,p)*leftSigma + 
							4.*spec->N(rCenter,p)*centerSigma +
							spec->N(rRight,p)*rightSigma);
						ASSERT_VALID_NO(sum);
					}
					m_sigmas(p)[nS][nratio]+=sum;
				}
			}
			centerSigma = m_difSigma->Sigma(S,ratio/ss2);
			rightSigma = leftSigma;
		}
	}
	for(nS = 0; nS<m_nn; nS++)
		for(nratio = nn-1; nratio>=0; nratio--)
		{
			FOR_ALL_REAL_PARTICLES_INVOLVED(particle)
			{
				if(!m_sigmas[particle])
					continue;
				m_sigmas(particle)[nS][nratio]*=integrationMultiplier;
			}
		}
	isPreLinkInitDone = 1;
	PostInit();
	isReady=1;
}

void CSecondaryDifSigma::AddSecondaryLink(TParticle _from, TParticle _to)
{
	ASSERT(isPreLinkInitDone);
	//links may be set only from PostInit() method
	//or after calling Init();

	m_sigmas[_from] = m_sigmas[_to];
}

void CSecondaryDifSigma::AddLink(TParticle _from, TParticle _to)
{
	ASSERT(isPreLinkInitDone);
	//links may be set only from PostInit() method
	//or after calling Init();

	m_sigmas[_from] = m_sigmas[_to];
	m_directRates[_from]=m_directRates[_to];
}

void CSecondaryDifSigma::MakeAntiparticlesLikeParticles()
{
	ASSERT(isPreLinkInitDone);
	//links may be set only from PostInit() method
	//or after calling Init();

	FOR_ALL_REAL_PARTICLES_INVOLVED(i)
	{
		TParticle from = (TParticle)i;
		TParticle to = ParticleData::Antiparticle2Particle(from);
		if(from!=to)
			AddLink(from,to);
	}
}

void CSecondaryDifSigma::DeleteSecondaryLink(TParticle _from)
{
	ASSERT(isPreLinkInitDone);
	//links may be set only from PostInit() method
	//or after calling Init();

	m_sigmas[_from]=NULL;
}

void CSecondaryDifSigma::MoveSecondaryLink(TParticle _from, TParticle _to)
{
	ASSERT(isPreLinkInitDone);
	//links may be set only from PostInit() method
	//or after calling Init();

	m_sigmas[_to]=m_sigmas[_from];
	m_sigmas[_from]=NULL;
}

double CSecondaryDifSigma::CalculateSigma_tot(int _s)
{
	const double Emin = Ranges().Emin();
	const double Emax = Ranges().Emax();
	const double ss1 = BC().ss1;
	double ratio = Emin/Emax*ss1;
	double multiplier = log(10.)/BC().s;
	double result = 0.;
	for(int r=0;r<m_nn;r++,ratio*=ss1)
		result+=(ratio*m_primarySigma[_s][r]);
	return result*multiplier;
}
