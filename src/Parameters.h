#ifndef IPAREADER_INCLUDED
#define IPAREADER_INCLUDED

#include "Addfunc.h"

#define DEFAULT_STR_PAR_MAX_LENGTH 1024

enum EParameterType{
	doubleE = 0,//order of elements in the enum should be synchronized with sTypeNames array values
	intE,
	boolE,
	switchE,
	stringE,

	endParameterTypeE
};

class IParReader
{
public:

	virtual ~IParReader(){}

	virtual bool parse()=0;

	virtual bool parExists(const char* theParName, EParameterType theParType) = 0;

	virtual bool readDoublePar( const char* aParName,
				   double& aPar)=0;

	virtual bool readIntPar( const char* aParName,
				   int& aPar)=0;

	virtual bool readBoolPar( const char* aParName,
				   bool &aPar)=0;

	virtual bool readSwitchPar( const char* aParName,
				   int& aPar)=0;
				   
	virtual bool readStringPar( const char* aParName,
				   char* aPar, int maxLength = DEFAULT_STR_PAR_MAX_LENGTH)=0;
	/*
	inline bool readboolPar( const char* aParName,
					   bool &aPar)
	{
		bool par = aPar;
		bool result = readboolPar(aParName, par);
		aPar = par != 0;
		return result;
	}*/

	//possible switch values expected:   aMinValue:aEndValue-1
	inline bool readSwitchPar(std::string aParName, void* aPar, int aEndValue, int aMinValue=0)
	{
		int* pPar = (int*)aPar;
		int par = *pPar;
		if(!readSwitchPar(aParName.c_str(), par))
			return false;
		if(par<aMinValue || par >= aEndValue)
			ThrowError("Unexpected parameter value " + aParName + "=" + ToString(par));
		*pPar = par;
		return true;
	}

	inline bool readStringPar( std::string aParName,
							   std::string& aPar, int maxLength = DEFAULT_STR_PAR_MAX_LENGTH){
		SafePtr<char> val = new char[maxLength];
		ASSERT(aPar.length()<maxLength);
		strncpy(val,aPar.c_str(),maxLength);
		bool result = readStringPar(aParName.c_str(),(char*)val, maxLength);
		if(result)
			aPar = val;
		return result;
	}
};

//base class for all classes using parameters read from config file
//parameters should be stored in static fields
//constructor is expected to read parameters only once
//This class has constructor which just checks that Parameters::fReader is initialized
class Parameters
{
public:
	static void SetReader(IParReader* aReader);
protected:
	Parameters();
	static inline IParReader* Reader() { return fReader; }
private:
	static IParReader* fReader;
};

#define READ_DOUBLE_SETTING(setting) Reader()->readDoublePar(#setting,setting)
#define READ_INT_SETTING(setting) Reader()->readIntPar(#setting,setting)
#define READ_BOOL_SETTING(setting) Reader()->readBoolPar(#setting,setting)
#define READ_STRING_SETTING(setting) Reader()->readStringPar(#setting,setting)
#define READ_SWITCH_SETTING(setting, maxValue) Reader()->readSwitchPar(#setting,&setting,maxValue)

class IParWriter
{
public:
	virtual bool print()=0;

	virtual bool writePar( const char* theParName,
				   double thePar, const char *theComment=0)=0;

	virtual bool writePar( const char* theParName,
				   int thePar, const char *theComment=0)=0;

	virtual bool writeBoolPar( const char* theParName,
				   bool thePar, const char *theComment=0)=0;

	virtual bool writeSwitchPar( const char* theParName,
				   int thePar, const char *theComment=0)=0;

	virtual	bool writePar( const char* theParName,
				   const char* thePar, const char *theComment=0)=0;

	virtual void WriteLine(const char* theComment = 0)=0;

	virtual bool SupportsPositioning() = 0;

	inline static IParWriter* Instance() { return fInstance; }
	inline static void Set(IParWriter* aValue){ fInstance = aValue; }
private:
	static IParWriter* fInstance;
};

#endif//end of file
