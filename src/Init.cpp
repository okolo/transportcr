// Init.cpp: implementation of the Initializer class.
//
//////////////////////////////////////////////////////////////////////

#include "Init.h"
#include "Addfunc.h"

Initializer globalInit;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
bool Initializer::isDone = false;

Initializer::Initializer()
{
	if (!isDone)
	{
		init();
		isDone = true;
	}
}

/**
Place here all program initialisation code
*/
void Initializer::init()
{
	initDirStructure();
}


void Initializer::initDirStructure()
{
	Mkdir(ALL_PLT_DIR);//this will create directory only if it doesn't exist yet
	Mkdir(PLT_DIR_LINK);
}
