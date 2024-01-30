// Init.h: interface for the Init class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INIT_H__767ECF0B_CC80_4C87_BC05_F69FFE597940__INCLUDED_)
#define AFX_INIT_H__767ECF0B_CC80_4C87_BC05_F69FFE597940__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Initializer  
{
public:
	Initializer();
private:
	static void initDirStructure();
	static void init();
// 
static bool isDone;
};

#endif // !defined(AFX_INIT_H__767ECF0B_CC80_4C87_BC05_F69FFE597940__INCLUDED_)
