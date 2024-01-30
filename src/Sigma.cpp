#include "Sigma.h"
#include "Addfunc.h"
#include "FilePtr.h"
#include "Units.h"

//print sigma (in sm^2) to file _n - number of points
void CSigma::Print(string _fileName, double _sMin, int _n, double _multipier) const
{
	double S = _sMin;
	double eMult = units.Eunit/units.outEunitMeV;
	eMult*=eMult;
	FILE* file = Fopen(_fileName,"wt");
	for(int i=0;i<_n;i++,S*=BC().ss1)
	{
		double sigmaVal=Sigma(S)*_multipier;
		sigmaVal*=units.sigmaUnit;
		fprintf(file,"%lg %lg\n",S*eMult,sigmaVal);
	}
	fclose(file);
}

double CSigma::MeanSigma(double sLeft, double sRight, int accuracy) const
{
	ASSERT(accuracy>0);
	double S = sLeft;
	double left = Sigma(S);
	double sum = 0.;
	double dS = (sRight - sLeft)/((double)accuracy);
	for(int i=0;i<accuracy;i++)
	{
		sum += 4.*Sigma(S+0.5*dS)+left;
		S+=dS;
		left = Sigma(S);
		sum += left;
	}
	
	ASSERT(sum>=0.);
	return sum/6./accuracy;
}

double CDifSigma::SigmaTotal(double _s, int _accuracy) const
{
	ASSERT(_accuracy>=2);
	double dr=1./((double)_accuracy);
	double r = dr;
	double leftVal = Sigma(_s,r);
	double result=0.;
	do
	{
		result+=leftVal;
		result+=4.*Sigma(_s,r+0.5*dr);
		r+=dr;
		leftVal = Sigma(_s,r);
		result+=leftVal;
	}while(r<1.);
	result*=(dr/6.);
	return result;
}

double CDifSigma::SigmaTotalE(double _s, int _accuracy, int noOfSecondaries) const
{
	ASSERT(_accuracy>=2);
	double dr=1./((double)_accuracy);
	double r = dr;
	double dr2 = 0.5*dr;
	double leftVal = r*Sigma(_s,r);
	double result=0.;
	do
	{
		result+=leftVal;
		r+=dr2;
		result+=4.*Sigma(_s,r)*r;
		r+=dr2;
		leftVal = Sigma(_s,r)*r;
		result+=leftVal;
	}while(r<1.);
	result*=(dr/6.*noOfSecondaries);
	return result;
}

void CDifSigma::printSigmaTot(string _fileName, double _sMin, int _n, int _accuracy) const
{
	CFilePtr file(Fopen(_fileName,"wt"));

	double S = _sMin;
	double eMult = units.Eunit/units.outEunitMeV;
	eMult*=eMult;
	for(int i=0;i<_n;i++,S*=BC().ss1)
	{
		double sigmaVal=SigmaTotal(S,_accuracy);
		sigmaVal*=units.sigmaUnit;
		fprintf(file,"%lg %lg\n",S*eMult,sigmaVal);
	}
}

void CDifSigma::printSigmaTotE(string _fileName, double _sMin, int _n, int _accuracy, int noOfSecondaries) const
{
	CFilePtr file(Fopen(_fileName,"wt"));

	double S = _sMin;
	double eMult = units.Eunit/units.outEunitMeV;
	eMult*=eMult;
	for(int i=0;i<_n;i++,S*=BC().ss1)
	{
		double sigmaVal=SigmaTotalE(S,_accuracy, noOfSecondaries);
		sigmaVal*=units.sigmaUnit;
		fprintf(file,"%lg %lg\n",S*eMult,sigmaVal);
	}	
}
