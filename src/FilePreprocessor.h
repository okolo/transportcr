// FilePreprocessor.h: interface for the CFilePreprocessor class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEPREPROCESSOR_H__12762CAD_70F8_11D4_A2CC_005004BE1D15__INCLUDED_)
#define AFX_FILEPREPROCESSOR_H__12762CAD_70F8_11D4_A2CC_005004BE1D15__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "stdio.h"

class CFilePreprocessor  // this class is used to remove comments from parameter files
{
protected:
	FILE* m_PreprocessedFile;
	char* m_PreprocessedFileName;
	char m_commentSimbol;
public:
	CFilePreprocessor(char commentSimbol='#');
	int Open(const char* fileName);
	void Close();
	virtual ~CFilePreprocessor();
	operator FILE*(){return m_PreprocessedFile;};

};

#endif // !defined(AFX_FILEPREPROCESSOR_H__12762CAD_70F8_11D4_A2CC_005004BE1D15__INCLUDED_)
