// TableReader.cpp: implementation of the CTableReader class.
//
//////////////////////////////////////////////////////////////////////

#include "TableReader.h"
#include "FilePtr.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTableReader::CTableReader(std::string aFilePath, int _noOfColumns):
m_columns(_noOfColumns),
fFilePath(aFilePath)
{
	CFilePtr file(FopenL(aFilePath,"rt"));
	Init(file, _noOfColumns);
}

CTableReader::CTableReader(FILE* _file, int _noOfColumns):
m_columns(_noOfColumns)
{
	Init(_file, _noOfColumns);
}

void CTableReader::Init(FILE* _file, int _noOfColumns)
{
	const double Nan = -1.3e300;
	int bufSize = _noOfColumns*64;
	char* lineBuffer = new char[bufSize];
	int i;
	for(i=0;i<_noOfColumns;i++)
		m_columns[i] = new CVarVector;
	//double value;

	std::string formatStr = "";
	for(i=0;i<_noOfColumns;i++)
		formatStr = formatStr + " %lg";

	for(int line = 1; fgets(lineBuffer,bufSize,_file); line++)//read line
	{
		if(lineBuffer[0]=='#')
			continue;//skip comments

		std::istringstream theLine(lineBuffer);

		for(i=0;i<_noOfColumns;i++)
		{
			double value = Nan;
			theLine >> value;
			if(value == Nan)
			{
				if(i>0)//empty strings are allowed
					ThrowError("CTableReader invalid '" + fFilePath + "' file format (line #" + ToString(line) + ")");
				break;
			}
			m_columns[i]->add(value);
		}
	}
}

CTableReader* CTableReader::ShiftCol(int aCol, double aShift){
    CVarVector& col = getColumn(aCol);
    for(int i=col.length()-1; i>=0; i--)
        col[i] += aShift;
    return this;
}

CTableReader* CTableReader::MultiplyCol(int aCol, double aMult){
    CVarVector& col = getColumn(aCol);
    for(int i=col.length()-1; i>=0; i--)
        col[i] *= aMult;
    return this;
}
