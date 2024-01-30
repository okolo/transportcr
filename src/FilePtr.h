// FilePtr.h: interface for the CFilePtr class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEPTR_H__9C824720_C230_11D5_885E_444553540001__INCLUDED_)
#define AFX_FILEPTR_H__9C824720_C230_11D5_885E_444553540001__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>

class CFilePtr  
{
public:
	CFilePtr(FILE* file = NULL):m_file(file){}
	operator FILE*(){return m_file;};
	virtual ~CFilePtr();
protected:
	FILE* m_file;
};

#endif // !defined(AFX_FILEPTR_H__9C824720_C230_11D5_885E_444553540001__INCLUDED_)
