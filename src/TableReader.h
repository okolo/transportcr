// TableReader.h: interface for the CTableReader class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(TABLEREADER_H__INCLUDED_)
#define TABLEREADER_H__INCLUDED_

#include "Vector.h"
#include <stdio.h>

class CTableReader  
{
public:
	CTableReader(std::string aFilePath, int _noOfColumns);
	CTableReader(FILE* _file, int _noOfColumns);
	inline CVarVector& getColumn(int i){return *(m_columns[i]);};
	inline const CVarVector& getColumn(int i) const {return *(m_columns[i]);};
	inline int ColumnCount() const { return m_columns.length(); }
    CTableReader* ShiftCol(int aCol, double aShift);
    CTableReader* MultiplyCol(int aCol, double aMult);
protected:
	void Init(FILE* _file, int _noOfColumns);
	CAutoDeletePtrArray<CVarVector> m_columns;
	std::string						fFilePath;
};

#endif // !defined(TABLEREADER_H__INCLUDED_)
