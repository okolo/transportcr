#include "Addfunc.h"
#include "ScanInfoWritter.h"

using namespace std;


ScanInfoWritter::~ScanInfoWritter(void)
{
	fOutput.flush();
}

template<class T> void PrintParams (std::vector<KeyValuePair<T> >& aVector, std::ostream& aOut, bool aPrintNamesInstead = false)
{
	//sortparams alphabetically before printing
	std::sort(aVector.begin(), aVector.end());
	for(unsigned int i=0; i<aVector.size(); i++) 
    {
		if(aPrintNamesInstead)
			aOut << aVector[i].key.c_str() <<" ";
		else
			aOut << aVector[i].value <<" ";
    }
	aOut << std::endl;
}

bool ScanInfoWritter::print()
{
	if(fFilteredParams.size() == 0)
	{
		PrintParams(fDoubleParams, fOutput, fPrintNamesInstead);
		PrintParams(fIntParams, fOutput, fPrintNamesInstead);
		PrintParams(fboolParams, fOutput, fPrintNamesInstead);
		PrintParams(fSwitchParams, fOutput, fPrintNamesInstead);
	}
	else
	{//print only given params
		int iMax = fFilteredParams.size();
		for(int i=0; i<iMax; i++)
		{
			if(i > 0)
				fOutput << "\t";
			fOutput << fFilteredParams[i];
		}
		fOutput.flush();
	}
	return true;
}

void ScanInfoWritter::setFilter(char *argv[], int aLength)
{
	fFilter.clear();
	fFilteredParams.clear();
	for(int i=0; i < aLength; i++)
	{
		fFilter.push_back(argv[i]);
		fFilteredParams.push_back("UNKNOWN");
	}
}

template<class T> void FilterPar(const char* name, const T& aPar, 
	std::vector <std::string >& aFilter, std::vector <std::string >&	aFilteredParams)
{
	for(int i=aFilter.size()-1; i>=0; i--)
	{
		if(aFilter[i] == name)
		{
			aFilteredParams[i] = ToString(aPar);
			break;
		}
	}
}

bool ScanInfoWritter::writePar( const char* theParName,
	   double thePar, const char *theComment)
{
	FilterPar(theParName, thePar, fFilter, fFilteredParams);
	KeyValuePair<double> pair(theParName, thePar);
	fDoubleParams.push_back(pair);
	return true;
}

bool ScanInfoWritter::writePar( const char* theParName,
	   int thePar, const char *theComment)
{
	FilterPar(theParName, thePar, fFilter, fFilteredParams);
	KeyValuePair<int> pair(theParName, thePar);
	fIntParams.push_back(pair);
	return true;
}

bool ScanInfoWritter::writeBoolPar( const char* theParName,
	   bool thePar, const char *theComment)
{
	FilterPar(theParName, thePar, fFilter, fFilteredParams);
	KeyValuePair<bool> pair(theParName, thePar!=0);
	fboolParams.push_back(pair);
	return true;
}

bool ScanInfoWritter::writeSwitchPar( const char* theParName,
	   int thePar, const char *theComment)
{
	FilterPar(theParName, thePar, fFilter, fFilteredParams);
	KeyValuePair<int> pair(theParName, thePar);
	fSwitchParams.push_back(pair);
	return true;
}

bool ScanInfoWritter::writePar( const char* theParName,
	   const char* thePar, const char *theComment)
{
	FilterPar(theParName, thePar, fFilter, fFilteredParams);
	return true;//don't save string params
}
	   
void ScanInfoWritter::WriteLine(const char* theComment)
{
	//do nothing
}
