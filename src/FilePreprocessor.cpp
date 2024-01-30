// FilePreprocessor.cpp: implementation of the CFilePreprocessor class.
//
//////////////////////////////////////////////////////////////////////

#include "FilePreprocessor.h"
#include <string.h>
#include <assert.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CFilePreprocessor::CFilePreprocessor(char commentSimbol)
{
	m_commentSimbol=commentSimbol;
	m_PreprocessedFile=NULL;
	m_PreprocessedFileName=NULL;
}

void CFilePreprocessor::Close()
{
	if(m_PreprocessedFile!=NULL)
	{
		fclose(m_PreprocessedFile);
		remove(m_PreprocessedFileName);
	};
	m_PreprocessedFile=NULL;
	delete m_PreprocessedFileName;
}

CFilePreprocessor::~CFilePreprocessor()
{
	Close();
}

int CFilePreprocessor::Open(const char* fileName)
{
	Close();
	FILE* originalFile=fopen(fileName,"rt");
	if(originalFile==NULL)
		return 0;

	m_PreprocessedFileName=new char[strlen(fileName)+15];
	strcpy(m_PreprocessedFileName,fileName);
	strcat(m_PreprocessedFileName,"_prepr_tmp_");

	FILE* preprocessedFile=fopen(m_PreprocessedFileName,"wt");
	assert(preprocessedFile!=NULL);
	int isComment=0;
	for(int ch=fgetc(originalFile);ch!=EOF;ch=fgetc(originalFile))
	{
	  //		char c=ch;
		if(isComment)
		{
			if(ch=='\n')
			{
				isComment=0;
				fputc(ch,preprocessedFile);
			};
		}
		else
		{
			if(ch==m_commentSimbol)
				isComment=1;
			else
				fputc(ch,preprocessedFile);
		}
	}
	fclose(originalFile);
	fclose(preprocessedFile);
	m_PreprocessedFile=fopen(m_PreprocessedFileName,"rt");
	if(m_PreprocessedFile==NULL)
		return 0;
	return 1;//true
}
