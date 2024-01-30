/*
 * Sarkar1005Background.cpp
 *
 *  Created on: May 19, 2011
 *      Author: ok
 */

#include "Sarkar1005Background.h"
#include "FilePtr.h"
#include "TableReader.h"

Sarkar1005Background::Sarkar1005Background(SarkarBackgrType aType):
f0(0),
fEvol(0)
{
	const char* fileName = aType==SarkarRadio ? "radioSarkar" : "iroSarkar";
	std::string dir = DATA_DIR "Sarkar1005Background" DIR_DELIMITER_STR;
	CFilePtr dataFilePtr(Fopen(dir + fileName, "rt"));
	f0 = new CLogScaleLinearFunc(new CTableReader(dataFilePtr, 2));
	CFilePtr dataFilePtrEvol(Fopen(dir + "evolutionSarkar", "rt"));
	CTableReader* reader = new CTableReader(dataFilePtrEvol, 2);
	fMaxZ = reader->getColumn(0)[reader->getColumn(0).length()-1];
	fEvol = new CLinearFunc(reader);
}

Sarkar1005Background::~Sarkar1005Background() {
	delete f0;
	delete fEvol;
}

double Sarkar1005Background::F(double E, double z)
{//fEvol->f(0) === 1
//taken from formula (C3) of http://lanl.arxiv.org/abs/0902.3993v2
//but divided by factor (1+z)^3 as required by definition of IBackgroundSpectrum::F
//and multiplied by (1+z) since the output is n(k)*k and not n(k)
	return pow(1.+z,-3)*fEvol->f(z)*f0->f(E/(1.+z));
}

double Sarkar1005Background::MaxE(double aZmax) const
{
	return f0->lastX()*(1.+MIN(fMaxZ, aZmax));
}

double Sarkar1005Background::MinE(double aZmax) const
{
	return f0->firstX();
}
