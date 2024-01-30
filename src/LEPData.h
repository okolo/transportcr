// LEPData.h: interface for the CLEPData class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LEPDATA_H__0EAC35C6_46DF_11D5_B6B3_C23BACEB308A__INCLUDED_)
#define AFX_LEPDATA_H__0EAC35C6_46DF_11D5_B6B3_C23BACEB308A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DataReader.h"

class CLEPData : public CTableWithHeaderReader  
{
public:
	void Save();
	CLEPData();
	virtual ~CLEPData();
	enum EParticles
	{
		electronsE = 0,
		photonsE,
		protonsE,
		neutronsE,
		neutrinosE,
		endE
	};
	double N(EParticles particle,double _x);

protected:
	virtual ECompleteStatus readHeaderLine(const char* theString);
	virtual bool readDataLine(const char* theString);
	virtual bool testData();
	virtual void processData();

private:
	double* x;
	double** y;
	unsigned int m_accuracy;
	unsigned int m_noOfDecays;
	unsigned int m_currentPoint;
	double m_controlSum;
};

#endif // !defined(AFX_LEPDATA_H__0EAC35C6_46DF_11D5_B6B3_C23BACEB308A__INCLUDED_)
