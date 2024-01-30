// DataReader.cpp: implementation of the CDataReader class.
//
//////////////////////////////////////////////////////////////////////

#include "DataReader.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTableWithHeaderReader::CTableWithHeaderReader(const char* theFileName, char theCommentSimbol):
m_commentSymbol(theCommentSimbol),
m_curHeaderLine(0),
m_curDataLine(0)
{
	m_file = fopen(theFileName,"rt");
}

CTableWithHeaderReader::~CTableWithHeaderReader()
{

}

bool CTableWithHeaderReader::readLine()
{
	for(int i=0;i<BUF_LENGTH;i++)
	{
		int ch = fgetc(m_file);
		if(ch == EOF)
		{
			if(m_buffer[0] == m_commentSymbol)
				i = 0;
			m_buffer[i] = '\0';
			return i;
		}
		if(ch == '\n')
		{
			if((i==0)||(m_buffer[0] == m_commentSymbol))
			{
				i = -1;
				continue;//skipping empty strings and comments
			}
			else
			{
				m_buffer[i] = '\0';
				return true;
			}
		}
		else
			m_buffer[i] = ch;
	}
	fprintf(stderr,"CDataReader : buffer overflow");
	Exit(retBufferOverflow);
	return false;
}

bool CTableWithHeaderReader::read()
{
	if(!m_file)
		return false;
	//	bool result = false;
	while(readLine())
	{
		m_curHeaderLine++;
		ECompleteStatus status = readHeaderLine(m_buffer);
		switch(status)
		{
			case notCompletedE: continue;
			case failedE: return false;
			case completedE:;
		}
		break;
	}
	while(readLine())
	{
		m_curDataLine++;
		if(!readDataLine(m_buffer))
			return false;
	}
	if(!testData())
		return false;
	processData();
	fclose(m_file);
	m_file = NULL;
	return true;
}

/// reads number from standard file header line "data length <number>"
/// here number must indicate number of data lines below
bool CTableWithHeaderReader::readDataLength(int& length, const char* theString)
{
	const char format[] = " data length %d";
	return (sscanf(theString,format,&length)==1);
}


/// parameter vector must not be created previously (using CVector::Create() method)
/// default delimiter string " /t/r/n"
bool CTableWithHeaderReader::readVector(CVector& v, const char* theString, const char* delimiterStr)
{
	bool result = false;
	char* stringCopy = 0;
	double* tmpArray = 0;
	try{
		const char defaultDelimStr[] = " \t\n\r";
		if (delimiterStr==NULL)
		{
			delimiterStr = defaultDelimStr;
		}

		int length = strlen(theString);
		stringCopy = new char[length + 1];
		strcpy(stringCopy, theString);

		tmpArray = new double[length/2 + 1]; // maximal possible number of numbers stored in the string

		int i=0;
		char* curStr = strtok(stringCopy,delimiterStr);

		// calculating number of entries
		while(curStr)
		{
			if (sscanf(curStr,"%lg",tmpArray+i)!=1)
				ThrowError("number format error");
			curStr = strtok(NULL, delimiterStr);
			i++;
		}
		if(i==0)
		{
			ThrowError("no tokens");
		}

		v.create(i);
		v.copy(tmpArray);

		result = true;
	}
	catch(const char* aError)
	{
		//error saved and can be retrieved using LastError() method
	}
	delete[] stringCopy;
	delete[] tmpArray;
	return result;
}

double CTableWithHeaderReader::yValue(double xValue, double *xArray, double *yArray, int xArraySize, double leftValue, double rightValue)
//find the nearest to xValue x-array element (increasing x-array supposed)
//like CVector::findX function, and than using 1-st order approximation to find y value
{
	int left = 0;
	int right = xArraySize-1;

	if(xValue>xArray[right])
		return rightValue;
	if(xValue<xArray[left])
		return leftValue;

	for(int i=(left+right)/2;right-left>1;i=(left+right)/2)
	{
		if(xValue>xArray[i])
			left=i;
		else
			right=i;
	}//finding nearest point

	ASSERT((right - left) == 1);

	double result = yArray[left]+(xValue-xArray[left])*(yArray[right]-yArray[left])/(xArray[right]-xArray[left]);
	return result;
}

void CTableWithHeaderReader::inverseArray(double *pArray, unsigned int size)
{
	unsigned int maxIndex = size/2;
	unsigned int i1,i2;
	for(i1=0,i2=size-1;i1<maxIndex;i1++,i2--)
	{
		double mem = pArray[i1];
		pArray[i1] = pArray[i2];
		pArray[i2] = mem;
	}
}
