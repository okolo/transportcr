#include <string>
#include <iostream>
#include <iomanip>
#include "Xmax.h"
#include "TableReader.h"
#include "TableFunc.h"
#include "TParticle.h"
#include "Vector.h"
#include "FilePreprocessor.h"
#include "VersionInfo.h"

/// propagation stub (to avoid link problems)
std::string plt_local_dir(PLT_DIR_LINK);
using namespace std;

static void printXmaxUsage()
{
	VersionInfo::PrintVersionInfo(cerr, "Xmax");
	cerr << "Xmax: utility to calculate Xmax and standard deviation of Xmax based on the spectrum and hadron model data\n" <<
		"usage Xmax <XmaxDataDir> Emin Emax <spectrum_file>\n" <<
		"example: ./Xmax tables/Xmax/QGSJET01 1e19 1e20 results/myRun/uniform/all\n\n" <<
		"Xmax --printPatch file.diff\n"
		<< "print git patch file (useful for dirty version builds)" << endl;
}

int main(int argc, char* argv[])
{
	if(argc==3 && strcmp(argv[1],"--printPatch")==0)
	{
		VersionInfo::PrintDiffFile(argv[2]);
		return 0;
	}
	if(argc != 5)
	{
		printXmaxUsage();
		return 1;
	}
	double Emin, Emax;
	if(sscanf(argv[2], "%lg", &Emin)!=1)
	{
		cerr << "invalid Emin parameter (not a number)\n";
		printXmaxUsage();
		return 1;
	}
	if(sscanf(argv[3], "%lg", &Emax)!=1)
	{
		cerr << "invalid Emax parameter (not a number)\n";
		printXmaxUsage();
		return 1;
	}
	if (Emax < Emin)
	{
		cerr << "invalid (Emin, Emax) range: (Emax < Emin)\n";
		printXmaxUsage();
		return 1;	
	}
	try
	{
		XMaxCalc calc(argv[1], Emin, Emax, argv[4]);
		calc.run();
	}
	catch(int aRetCode)
	{
		return aRetCode;
	}
	catch(string aErrorMessage)
	{
		cerr << aErrorMessage << endl;
		return 2;
	}
	catch(const char* aErrorMessage)
	{
		cerr << aErrorMessage << endl;
		return 2;
	}
	return 0;
}

//optimal scale functions below are used to achieve better accuracy
//using CLinearFunc class
//optimal scale is choosen so that the quantity of intereset
//depends approximately linearly on the argument

inline double optimalScaleForXmax(double A)
{//unger2.pdf (slide 8): <Xmax> ~ lg(E/A) + const
	return log10(A);
}

inline double optimalScaleForSigma(double A)
{//NOTE!!!: it is supposed farther in programm that the optimal A scale for sigma is decreasing function of A

	return pow(A, -0.2);//looking at fig. on slide 9 of unger2.pdf
	//(QGSJET01 points for p,He,N and Fe at E=1e17 eV are 70,50,36,24 which
	//is very well fitted by A^(-1/5) power low
}

XMaxCalc::XMaxCalc(string aTableDir, double aEmin, double aEmax, string aSpectrumFile):
fXmax(57),//fXmax[0] will not be used
fXmaxSigma(57),//fXmaxSigma[0] will not be used
fSpectrum(0),
fFractions(57),
fCurXmax(0),
fCurXmaxSigma(0),
fEmin(aEmin),
fEmax(aEmax)
{
	//reading spectrum
	{
		CFilePreprocessor specFile;//removing header staring with '#'
		if(!specFile.Open(aSpectrumFile.c_str()))
			throw "failed to open " + aSpectrumFile;

		fSpectrum = new CTableReader(specFile, EEndAllParticles+2);//Energy, spectra[EEndAllParticles], sum
		CVarVector& E = fSpectrum->getColumn(0);
		if(E[0]>aEmin)
			throw "spectrum Emin is higher than " + ToString(aEmin);
		//if(E[E.length()-1]<aEmax)
		//throw "spectrum Emax is lower than " + ToString(aEmax);
	}
	//reading tables
	aTableDir += DIR_DELIMITER_STR;
	CVarVector xPoints;
	CVarVector sPoints;
	for(int A=1; A<=56; A++)
	{
		string errorMessageStart = "File '";
		string errorMessageEnd = ToString(A) + "' not found in " + aTableDir;
		CFilePtr xFile = fopen((aTableDir + "x" + ToString(A)).c_str(), "rt");
		if(xFile!=NULL)
		{
			fXmax[A] = new CLinearFunc(new CTableReader(xFile, 2));//two column file expected
			fXmax[A]->SetExtension(ExtConst);
			xPoints.add(optimalScaleForXmax(A));
		}
		else if(A==1 || A==56)
			throw errorMessageStart + "x" + errorMessageEnd;//x1 and x56 should always be present

		CFilePtr sFile = fopen((aTableDir + "s" + ToString(A)).c_str(), "rt");
		if(sFile!=NULL)
		{
			fXmaxSigma[A] = new CLinearFunc(new CTableReader(sFile, 2));//two column file expected
			fXmaxSigma[A]->SetExtension(ExtConst);
			sPoints.add(optimalScaleForSigma(A));
		}
		//else if(A==1 || A==56)
		//	throw errorMessageStart + "x" + errorMessageEnd;//s1 and s56 should always be present
	}
	fXmaxA.create(xPoints.length());
	fXmaxValue.create(xPoints.length());
	fXmaxA = xPoints;
	if(sPoints.length()>=2)
	{
		fSigmaA.create(sPoints.length());
		fSigmaValue.create(sPoints.length());
		fSigmaA = sPoints;
		fSigmaA.invert();//optimal A scale for sigma is decreasing function of A
	}//else don't calculate sigma
}

XMaxCalc::~XMaxCalc()
{
	delete fSpectrum;
	delete fCurXmax;
	delete fCurXmaxSigma;
}

void XMaxCalc::run()
{
	CVarVector& E = fSpectrum->getColumn(0);
	int i;
	int iMax = E.length();
	for(i=0; E[i]<fEmin/1.001; i++);
	cout << setprecision (15);
	for(;i<iMax && E[i]<fEmax*1.001; i++)
	{
		//reading fractions
		double norm = 0.;
		for(int A=0; A<=56; A++)//including protons and neutrons
		{
			double flux = (fSpectrum->getColumn(ENeutron + A + 1))[i];//energy column is first
			fFractions[A] = flux;
			norm += flux;
		}
		if(norm==0.)//zero spectrum, don't print any data
			continue;
		norm = 1./norm;
		for(int A=0; A<=56; A++)//including protons and neutrons
		{
			fFractions[A] = fFractions[A]*norm;
		}
		fFractions[1] += fFractions[0];//from now on the protons and neutrons are counted together
		fFractions[0] = 0;
		//calculating Xmax and SigmaXmax
		SetCurrentE(E[i]);
		double xMax = 0;
		double xMaxSigma = 0;
		CalculateMeanXmaxAndSigmaXmax(xMax, xMaxSigma);

		//output
		cout << log10(E[i]) << "\t" << xMax << "\t" << xMaxSigma << endl;
	}
}

void XMaxCalc::CalculateMeanXmaxAndSigmaXmax(double& aXmaxOut, double& aSigmaXmaxOut)
{
	aXmaxOut = 0.;
	double RMS2 = 0.;//RMS^2 (RMS = root mean square)
	for(int A=1; A<=56; A++)
	{
		double f = fFractions[A];
		double x = Xmax(A);
		
		aXmaxOut += f * x;
		if(fCurXmaxSigma)//if sigma data is available
		{
			double s = XmaxSigma(A);
			double RMS2_A = s*s + x*x;
			RMS2 += f*RMS2_A;
		}
	}
	aSigmaXmaxOut = fCurXmaxSigma?sqrt(RMS2-aXmaxOut*aXmaxOut):0;
}

void XMaxCalc::SetCurrentE(double aE/*eV*/)
{
	double lgE = log10(aE);
	int iX = 0;
	int iS = fSigmaA.length()-1;
	for(int A=1; A<=56; A++)
	{
		if(fXmax[A]!=0)
		{
			fXmaxValue[iX] = fXmax[A]->f(lgE);
			iX++;
		}
		if(fXmaxSigma[A]!=0)
		{
			fSigmaValue[iS] = fXmaxSigma[A]->f(lgE);
			iS--;
		}
	}
	delete fCurXmax;
	fCurXmax = new CLinearFunc(fXmaxA, fXmaxValue);
	fCurXmax->SetExtension(ExtConst);
	if(fSigmaA.length()>=2)
	{
		delete fCurXmaxSigma;
		fCurXmaxSigma = new CLinearFunc(fSigmaA, fSigmaValue);
		fCurXmaxSigma->SetExtension(ExtConst);
	}
}

double XMaxCalc::Xmax(double aA)
{
	return fCurXmax->f(optimalScaleForXmax(aA));
}

double XMaxCalc::XmaxSigma(double aA)
{
	return fCurXmaxSigma->f(optimalScaleForSigma(aA));
}
