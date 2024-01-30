#ifndef ScanInfoWritter_INCLUDED
#define ScanInfoWritter_INCLUDED

#include "Parameters.h"
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>


template<class T>
class KeyValuePair
{
public:
	std::string key;
	T			value;
	KeyValuePair<T>(const char* aKey, T aValue)
	{
		key = aKey;
		value = aValue;
	}
	KeyValuePair<T>(const KeyValuePair<T>& aKeyValuePair)
	{
		key = aKeyValuePair.key;
		value = aKeyValuePair.value;
	}
	bool operator < (const KeyValuePair& _Right) const
	{
		return key.compare(_Right.key) < 0;
	}
};

class ScanInfoWritter :
	public IParWriter
{
public:
	ScanInfoWritter(std::ostream& aOutput, bool aPrintNamesInstead = false) : fOutput(aOutput), fPrintNamesInstead(aPrintNamesInstead){};

	virtual ~ScanInfoWritter(void);

	virtual bool print();

	bool printParamNames();

	virtual bool writePar( const char* theParName,
				   double thePar, const char *theComment=0);

	virtual bool writePar( const char* theParName,
				   int thePar, const char *theComment=0);

	virtual bool writeBoolPar( const char* theParName,
				   bool thePar, const char *theComment=0);

	virtual bool writeSwitchPar( const char* theParName,
				   int thePar, const char *theComment=0);
	
	virtual	bool writePar( const char* theParName,
				   const char* thePar, const char *theComment=0);
				   
	virtual void WriteLine(const char* theComment = 0);

	virtual bool SupportsPositioning(){return true;};

	void setFilter(char *argv[], int aLength);
private:
	std::ostream&						fOutput;
	bool								fPrintNamesInstead;
	std::vector <KeyValuePair<double> >	fDoubleParams;
	std::vector <KeyValuePair<int> >	fIntParams;
	std::vector <KeyValuePair<int> >	fSwitchParams;
	std::vector <KeyValuePair<bool> >	fboolParams;
	std::vector <std::string >			fFilter;
	std::vector <std::string >			fFilteredParams;
};


#endif //#ifndef ScanInfoWritter_INCLUDED
