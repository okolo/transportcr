#ifndef PRIMACK_IROSPECTRUM_ALREADY_INCLUDED
#define PRIMACK_IROSPECTRUM_ALREADY_INCLUDED


#include "TableBackground.h"
#include "Addfunc.h"
#include "Vector.h"

#define MAX_IRO_Z_ACCURASSY 6//changing this number isn't enough
//to change the accuracy (see readDataLine() and parseZvaluesString())

class CPrimackIROSpectrum : public CTableWithHeaderBackgr
{
public:
	double F(double E, double z);

	CPrimackIROSpectrum();
	virtual ~CPrimackIROSpectrum();
	double MinE(double aZmax) const;
	double MaxE(double aZmax) const;
	double MaxZ() const;
	
protected:

// from CDataReader
	virtual ECompleteStatus readHeaderLine(const char* theString);
	virtual bool readDataLine(const char* theString);
	virtual bool testData();
	virtual void processData();

	double F(double E, UINT iZ);// may throw "invalid argument"
// own functions
	double F(double E, UINT iZ, int* iE);	

private:
	bool parseZvaluesString(const char* theString);
	int m_zAc;
	int	m_eAc;
	CMatrix m_eArray;
	CMatrix m_nArray;
	CMatrix m_eAlphaArray;//energy derivative of the spectrum in logscale (or power law parameter in linear scale)
	CVector m_zArray;
	CVector m_eBuffer;
	CVector m_nBuffer;
};


#endif //#ifndef PRIMACK_IROSPECTRUM_ALREADY_INCLUDED
