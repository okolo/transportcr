#if !defined(SAFE_OUTPUT_H_INCLUDED)
#define SAFE_OUTPUT_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

typedef basic_ostream<char, char_traits<char> >& char_ostream_ref;

class CSafeOutput
{
public:

	CSafeOutput(){}

	CSafeOutput(string aFileName):
	  fOutput(aFileName.c_str())
	{
		
	}

	~CSafeOutput()
	{
		close();
	}

	inline void open(string aFileName)
	{
		fOutput.open(aFileName.c_str());
	}

	inline void close()
	{
		if(fOutput.is_open())
			fOutput.close();
	}

	inline void flush()
	{
		if(fOutput.is_open())
			fOutput.flush();
	}

	inline ofstream& Output()
	{
		return fOutput;
	}

private:
	ofstream fOutput;
};

	inline CSafeOutput& operator << (CSafeOutput& aOut, double aValue)
	{
		if(aOut.Output().is_open())
			aOut.Output() << aValue;
		return aOut;
	}

	inline CSafeOutput& operator << (CSafeOutput& aOut, int aValue)
	{
		if(aOut.Output().is_open())
			aOut.Output() << aValue;
		return aOut;
	}

	inline CSafeOutput& operator << (CSafeOutput& aOut, string aValue)
	{
		if(aOut.Output().is_open())
			aOut.Output() << aValue;
		return aOut;
	}

	inline CSafeOutput& operator << (CSafeOutput& aOut, const char* aValue)
	{
		if(aOut.Output().is_open())
			aOut.Output() << aValue;
		return aOut;
	}

	inline CSafeOutput& operator << (CSafeOutput& aOut, char_ostream_ref func(char_ostream_ref))
	{
		if(aOut.Output().is_open())
			aOut.Output() << func;
		return aOut;
	}

#endif // #if !defined(SAFE_OUTPUT_H_INCLUDED)