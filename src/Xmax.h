#ifndef XMAX_H_ALREADY_INCLUDED
#define XMAX_H_ALREADY_INCLUDED

#include <string>
#include "Vector.h"
#include "TableFunc.h"

class XMaxCalc
{
public:
	XMaxCalc(std::string aTableDir, double aEmin, double aEmax, std::string aSpectrumFile);
	~XMaxCalc();
	void run();
private:
	void CalculateMeanXmaxAndSigmaXmax(double& aXmaxOut, double& aSigmaXmaxOut);
	double Xmax(double aA);
	double XmaxSigma(double aA);
	void SetCurrentE(double aE/*eV*/);
private:
	CAutoDeletePtrArray<CLinearFunc>	fXmax;
	CAutoDeletePtrArray<CLinearFunc>	fXmaxSigma;
	CTableReader*						fSpectrum;
	CVector								fFractions;
	CLinearFunc*						fCurXmax;
	CLinearFunc*						fCurXmaxSigma;
	CVector								fXmaxA;
	CVector								fXmaxValue;
	CVector								fSigmaA;
	CVector								fSigmaValue;
	double								fEmin;
	double								fEmax;
};

#endif //#ifndef XMAX_H_ALREADY_INCLUDED