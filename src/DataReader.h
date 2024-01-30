// CTableWithHeaderReader.h: interface for the CTableWithHeaderReader class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATAREADER_H__D6683D61_4629_11D5_B6B2_93EA015A338A__INCLUDED_)
#define AFX_DATAREADER_H__D6683D61_4629_11D5_B6B2_93EA015A338A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Addfunc.h"
#include "Vector.h"

#define BUF_LENGTH 65536

class CTableWithHeaderReader  
{
public:
	static void inverseArray(double* pArray, unsigned int size);
	static double yValue(double xValue, double* xArray, double* yArray,int xArraySize, double leftValue=0, double rightValue=0);
	enum ECompleteStatus
	{
		notCompletedE = -1,
		failedE = 0,
		completedE = 1
	};
	CTableWithHeaderReader(const char* theFileName, char theCommentSimbol = '#');
	virtual bool read();
	virtual ~CTableWithHeaderReader();

	/// reads number from standard file header line "data length <number>"
	/// here number must indicate number of data lines below
	/// usually called from readHeaderLine
	/// parameter theString is header string passed to readHeaderLine
	static bool readDataLength(int& length, const char* theString);

	/// parameter vector must not be created previously (using CVector::Create() method)
	/// default delimiter string " /t/r/n"
	static bool readVector(CVector& vector, const char* theString, const char* delimiterStr=NULL);

protected:
	virtual ECompleteStatus readHeaderLine(const char* theString){return completedE;};

	//readDataLine returns false in case of error
	virtual bool readDataLine(const char* theString) = 0;
	virtual bool testData(){return true;};
	virtual void processData(){};

	char m_commentSymbol;
	int m_curHeaderLine;//starting from 1
	int m_curDataLine;//starting from 1
private:
	bool readLine();//returns false if EOF was reached

	char m_buffer[BUF_LENGTH];
	FILE* m_file;
};

#endif // !defined(AFX_DATAREADER_H__D6683D61_4629_11D5_B6B2_93EA015A338A__INCLUDED_)
