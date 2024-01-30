// Sigma.h: 
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SIGMA_H__6DF110E3_4AEB_11D5_B6B8_E8293AD9338A__INCLUDED_)
#define AFX_SIGMA_H__6DF110E3_4AEB_11D5_B6B8_E8293AD9338A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>

using namespace std;

#include "ParticleData.h"
#include "Ranges.h"


class ISigma
{
public:
	virtual double Sigma(double _s) const = 0;
};

class CSigma : public ISigma
{
public:
	double MeanSigma(double sLeft, double sRight, int accuracy) const;
	//print sigma (in sm^2) to file _n - number of points
	void Print(string _fileName, double _sMin, int _n, double _multipier=1.) const;
};

class CDifSigma
{
public:
	void printSigmaTot(string _fileName, double _sMin, int _n, int _accuracy) const;
	void printSigmaTotE(string _fileName, double _sMin, int _n, int _accuracy, int noOfSecondaries = 2) const;
	double SigmaTotal(double _s, int _accuracy) const;
	double SigmaTotalE(double _s, int _accuracy, int noOfSecondaries = 2) const;
	double Sigma(double _s,double _E, double _nE) const {return Sigma(_s,_nE/_E)/_E;};
	virtual double Sigma(double _s,double _ratio) const = 0;//sigma multiplied by E (initial particle energy)
};

class CDiscretDifSigma
{
public:
	CDiscretDifSigma(double _Smin, int _nn=Ranges().nE()):
	  m_Smin(_Smin),m_nn(_nn){m_ratioMin = Ranges().Emin()/Ranges().Emax()*BC().ss1;};
/**
*	here double s = s_min*BC().ss1^_s
*	double ratio = Emin/Emax*BC().ss1^(1+_ratio)
*/
	virtual double Sigma(int _s,int _ratio,TParticle _particle) const = 0;//sigma multiplied by E (initial particle energy)
protected:
	double m_Smin;
	int m_nn;//number of points on s scale
	double m_ratioMin;
};

#endif
