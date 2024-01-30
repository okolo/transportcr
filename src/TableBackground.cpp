#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <math.h>
#include "TableBackground.h"
#include "TableReader.h"
#include "FilePtr.h"
#include "Addfunc.h"

TableBackground::TableBackground(string aDir, const char** aFileList, bool aIsLogscaleY):
fMinE(-1.),
fMaxE(-1.),
fRangesZmax(-1)
{//the files should be ordered increasingly
	double lastZ = -1;
	CVarVector zVals;
	
	for(int i=0; aFileList[i]; i++)
	{
		double z = -1.;
		z = atof(aFileList[i]);
		if(z<0)
			throw "Invalid table background file name";
		if(z<lastZ)
			throw "Invalid table background file order";

                
		CFilePtr dataFilePtr(Fopen(DATA_DIR + aDir + DIR_DELIMITER_STR + aFileList[i], "rt"));
                CTableReader* reader = new CTableReader(dataFilePtr, 2);
                CVarVector& y = reader->getColumn(1);
                double leftValue = aIsLogscaleY?(y[0]-1000):0;
		double rightValue = aIsLogscaleY?(y[y.length()-1]-1000):0;
		CLinearFunc* f = new CLinearFunc(reader, leftValue, rightValue);
		fData.add(f);
		zVals.add(z);
	}
	fZ.create(zVals.length());
	fZ.copy(zVals);
	fMaxDataZ = fZ[fZ.length()-1];
}

/* spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
(must be multiplied by (1+z)^3 before substituting farther)*/
double TableBackground::F(double aE, double aZ)
{
	int iZ = fZ.findLeftX(aZ);
	if(iZ<0)
		return 0;
	double x = scaleX(aE, aZ);
	double leftY = fData(iZ)(x);
	double rightY = fData(iZ+1)(x);
	double zL = fZ[iZ];
	double zR = fZ[iZ+1];
	double y = leftY + (rightY-leftY)/(zR-zL)*(aZ-zL);
	return unscaleY(y, aE, aZ);
}

void TableBackground::InitRanges(double aZmax)
{
	fMinE = 1e300;
	int iZ = fData.length()-1;
	for(; iZ>0 && fZ[iZ-1]>=aZmax; iZ--);//find the range of z bins used

	for(; iZ>=0; iZ--)
	{
		double z = fZ[iZ];
		double enrgy = unscaleX(fData[iZ]->firstX(), z);
		if(enrgy<fMinE)
			fMinE = enrgy;
		if(enrgy>fMaxE)
			fMaxE = enrgy;
		enrgy = unscaleX(fData[iZ]->lastX(), z);
		if(enrgy<fMinE)
			fMinE = enrgy;
		if(enrgy>fMaxE)
			fMaxE = enrgy;
	}
	fRangesZmax = aZmax;
}

double TableBackground::MinE(double aZmax) const
{
	if(fRangesZmax!=aZmax)
		((TableBackground*)(this))->InitRanges(aZmax);
	return fMinE;
}

double TableBackground::MaxE(double aZmax) const
{
	if(fRangesZmax!=aZmax)
		((TableBackground*)(this))->InitRanges(aZmax);
	return fMaxE;
}

MatrixBackground::MatrixBackground(string aTableFile, bool aIsComoving, bool aExtendToZero):
		fZtoIndex(0),
		fLogEtoIndex(0),
		fExtendToZero(aExtendToZero),
		fIsComoving(aIsComoving),
		fBuffer(0),
		fMaxE(-1.),
		fMinE(-1.)
{
	const char* ch = aTableFile.c_str();
	const size_t bufSize = 1048576;//1M
	fBuffer = new char[bufSize];
	std::istream_iterator<double> eos;
	std::ifstream file(aTableFile.c_str());

	if(file.eof())
		ThrowError("Invalid background table file " + aTableFile);
	file.getline(fBuffer, bufSize);
	std::istringstream ist(fBuffer);
	std::istream_iterator<double> iit (ist);
	CVarVector Z,Zindex,LogE,LogEindex;
	for(int index=0; iit!=eos; iit++, index++)
	{
		Z.add(*iit);
		Zindex.add(index);
	}
	fZ.create(Z);
	fZindex.create(Zindex);
	double Epriv = 0;
	int line = 1;
	while(!file.eof())
	{
		line++;
		file.getline(fBuffer, bufSize);
		std::istringstream ist(fBuffer);
		std::istream_iterator<double> iit (ist);
		if(iit==eos)
			break;
		double E = *iit;
		if(E<=0 || E<Epriv)
			ThrowError("Invalid background table file " + aTableFile + ": invalid energy value in line " +
				ToString(line));
		Epriv = E;
		LogE.add(log(E));
		LogEindex.add(LogE.length()-1);
		CVector* pN = new CVector(fZ.length());
		fLogEdN_dE.add(pN);
		int recNo = 0;
		for(iit++; iit!=eos; iit++,recNo++)
		{
			if(recNo>=pN->length())
				continue;
			double dN_dE = *iit;
			if(dN_dE<0)
				ThrowError("Invalid background table file " + aTableFile + ": line " +
									ToString(line) + " contains " + ToString(recNo+1) +
									" negative value");
			(*pN)[recNo] = (dN_dE>0?log(E*dN_dE):-700.);
		}
		if(pN->length()!=recNo)
			ThrowError("Invalid background table file " + aTableFile + ": line " +
					ToString(line) + " contains " + ToString(recNo+1) +
					" records while " + ToString(pN->length()+1) + " expected");
	}
	fLogE.create(LogE);
	fLogEindex.create(LogEindex);
	fZtoIndex = new CDefaultTableFunc(fZ, fZindex, aExtendToZero ? 0. : -1., -1.);
	fLogEtoIndex = new CDefaultTableFunc(fLogE, fLogEindex, -1., -1.);
	fMinE = exp(fLogE[0]);
	fMaxE = exp(fLogE[fLogE.length()-1]);
	delete fBuffer;
	fBuffer = 0;
}

double MatrixBackground::MaxZ() const
{
	return fZ[fZ.length()-1];
}

double MatrixBackground::MaxE(double aZmax) const
{//todo: take into account aZmax
	return fMaxE;
}

double MatrixBackground::MinE(double aZmax) const
{//todo: take into account aZmax
	return fMinE;
}

MatrixBackground::~MatrixBackground()
{
	delete fBuffer;
	delete fZtoIndex;
	delete fLogEtoIndex;
}

/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
  (must be multiplied by (1+z)^3 before substituting farther)*/
double MatrixBackground::F(double aE, double aZ)
{
	double logE = log(aE);
	double indexE = fLogEtoIndex->f(logE);
	if(indexE<0.)
		return 0.;
	double indexZ = fZtoIndex->f(aZ);
	if(indexZ<0.)
		return 0.;
	int iE = indexE;
	double fracE = indexE-iE;
	int iZ = indexZ;
	double fracZ = indexZ-iZ;
	double result = (*fLogEdN_dE[iE])[iZ]*(1.-fracZ)*(1.-fracE);
	if(fracZ>0)
		result += ((*fLogEdN_dE[iE])[iZ+1]*fracZ*(1.-fracE));
	if(fracE>0)
	{
		result += ((*fLogEdN_dE[iE+1])[iZ]*(1.-fracZ)*fracE);
		if(fracZ>0)
			result += ((*fLogEdN_dE[iE+1])[iZ+1]*fracZ*fracE);
	}
	result = exp(result);
	if(!fIsComoving)
	{
		double z1=1.+aZ;
		result /= (z1*z1*z1);
	}
	return result;
}

double PlainTableBackground::F(double E, double z) {
	double mult = 1.;
	if(!fConstComovingDensity) {
        // since F(E,z) should return comoving concentration (i.e. result will be multiplied by (1+z)^3 later to obtain
        // physical concentration here we divide by (1+z)^3 to make physical concentration constant
		mult = 1/(1.+z);
		mult *= (mult*mult);
	}
    if(fEvolution){
        double F0 = fEvolution->F(E, 0);
        if(F0<=0.)
            ThrowError("Invalid custom background evolution: F(E=" + ToString(E) + " eV, z=0)=0");
        mult *= fEvolution->F(E, z)/F0;
    }
	return mult*fTable->f(E);
}

double PlainTableBackground::MaxZ() const {
	return fZmax;
}

double PlainTableBackground::MaxE(double aZmax) const {
	double maxE = fTable->lastX();
    if(fEvolution)
    {
        double maxEz = fEvolution->MaxE(aZmax);
        if(maxEz<maxE)
            return maxEz;
    }
    return maxE;
}

double PlainTableBackground::MinE(double aZmax) const {
	double minE = fTable->firstX();
    if(fEvolution)
    {
        double minEz = fEvolution->MinE(aZmax);
        if(minEz>minE)
            return minEz;
    }
    return minE;
}

PlainTableBackground::PlainTableBackground(std::string aTableFile, bool aConstComovingDensity, IBackgroundSpectrum* aEvolution, double aZmax) :
		fConstComovingDensity(aConstComovingDensity),
		fZmax(aZmax),
        fEvolution(aEvolution)
{
	fTable = new CLogScaleLinearFunc(new CTableReader(aTableFile,2));
}

bool PlainTableBackground::init()
{
    if(fEvolution){
        if(fEvolution->MaxE(0) < fTable->lastX() || fEvolution->MinE(0) > fTable->firstX())
            ThrowError("Invalid custom background evolution energy range");
        double evMaxZ = fEvolution->MaxZ();
        if(evMaxZ<fZmax)
            fZmax = evMaxZ;
    }
    return true;
}