#include "Dpp.h"
#include "TableReader.h"
#include "FilePtr.h"
#include "Ranges.h"
#include "Medium.h"
#include "Units.h"

namespace couplings{

void Dpp::Channel_gamma_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		aCoef.Add(iPrim, iPrim, -fCoupling.aDPP(iPrim));
}

void Dpp::Channel_gamma_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim, iSec;
	for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, 0.5*fCoupling.bDPP(iPrim,iPrim));
		for(iSec = 0; iSec<iPrim; iSec++)
			aCoef.Add(iSec, iPrim, fCoupling.bDPP(iSec,iPrim));
	}
}

Dpp::Dpp():
m_fi(0),
m_SigmaTot(0),
m_SigmaBackgroundIntegral(Ranges().nE()),
m_midE(Ranges().midE())
{
	AddChannel(new Channel_gamma_e(this, EElectron));
	AddChannel(new Channel_gamma_e(this, JOIN_ELECTRON_POSITRON ? EElectron : EPositron));
	AddChannel(new Channel_gamma_gamma(this));

	CFilePtr file(Fopen(DATA_DIR "DPP.tab","rt"));
	//tests for other tables
	//CFilePtr file(Fopen(DATA_DIR "dpp_0_126.txt","rt"));
	//CFilePtr file(Fopen(DATA_DIR "dpp_0_7924.txt","rt"));
	CTableReader* tr = new CTableReader(file, 2);
	CVarVector& xCol = tr->getColumn(0);
	CVarVector& yCol = tr->getColumn(1);
	
	int i=xCol.length();
	//if(i!=22)
	//	throw "invalid 'DPP.tab' file format (22 lines expected)";
	CLinearFunc f(tr);

	//normalizing
	double norm = f.integrate(xFunc);//energy
	//double norm = f.integrate();//probability
	if(norm <=0)
		throw "invalid data in 'DPP.tab'";
	norm = 1./norm;//energy
	//norm = 2./norm;//probability

	m_x.create(i);
	m_y.create(i);

	for(i--; i>=0; i--)
	{
		m_x[i] = xCol[i];
		m_y[i] = norm * yCol[i];
	}
	m_fi = new CLinearFunc(m_x, m_y);
	//norm = m_fi->integrate(xFunc);
	//norm = m_fi->integrate();

	m_sF = BC().sF;
	m_threshold = 4.*Me2;//here m_threshold == Sthreshold/4 (by definition), where S_threshold=16Me^2 CMF energy threshold for DPP
	InitSigmaTot();
}

void Dpp::InitSigmaTot()
{
	const double SigmaDPP=6.45e-30/units.sigmaUnit;
	CFilePtr file(Fopen(DATA_DIR "sigmaDpp.tab","rt"));
	CTableReader* tr = new CTableReader(file, 2);
	CVarVector& xCol = tr->getColumn(0);
	CVarVector& yCol = tr->getColumn(1);
	
	int i=xCol.length()-1;
	if(yCol[i]<=0)
		throw "invalid data in 'sigmaDpp.tab'";
	double norm = SigmaDPP/yCol[i];//normalizing last bin on asymptotic value
	for(; i>=0; i--)
	{
		xCol[i] *= Me2;//file column is in units of Me^2
		yCol[i] *= norm;//normalizing last bin on asymptotic value
	}

	m_SigmaTot = new CLinearFunc(tr, 0., SigmaDPP);
}

Dpp::~Dpp(void)
{
	delete m_fi;
	delete m_SigmaTot;
}

void Dpp::SetBackgrounds(const Medium& aPropagCoef)
{
	Coupling::SetBackgrounds(aPropagCoef);
	SetBackground(aPropagCoef.BackgroundIntegral());
}

void Dpp::SetBackground(const CBackgroundIntegral& aBackgr)
{
	m_Backgr = &aBackgr;
	double kMax = aBackgr.Kmax();
	int i = 0;
	for(i=m_midE.length()-1; i>=0; i--)
	{
		m_fParamE = m_midE[i];
		double kMin = m_threshold/m_fParamE;
		if(kMin>=kMax)
			break;
		m_SigmaBackgroundIntegral[i] = FuncUtils::integrateLogscale(kMin, kMax, this, 1000);
	}
	for(; i>=0; i--)
	{
		m_SigmaBackgroundIntegral[i] = 0.;
	}
}

double Dpp::f(double aX) const
{
	ASSERT(m_Backgr);
	return aX * m_Backgr->integral(aX) * m_SigmaTot->f(4.*aX*m_fParamE);
}

}//namespace couplings
