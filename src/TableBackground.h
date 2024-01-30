#if !defined(TableBackground_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
#define TableBackground_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_

#include "Background.h"
#include "Vector.h"
#include "TableFunc.h"
#include <string>
#include <vector>
#include "DataReader.h"
using namespace std;

class TableBackground :
	public IBackgroundSpectrum
{
public:
	//file names should be convertable to double (represent value of z)
	//last element in the array should be NULL
	//the files should be ordered increasingly
	TableBackground(string aDir, const char** aFileList, bool aIsLogscaleY);
	
	/* spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	(must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z);

	virtual double MaxZ() const { return fMaxDataZ; }

	//returns maximal background energy in eV
	double MaxE(double aZmax) const;

	//returns minimal background energy in eV
	double MinE(double aZmax) const;

protected:
	//convert aE to internal x scale
	virtual double scaleX(double aE/*eV*/, double aZ) = 0;

	//convert internal x scale to energy in eV
	virtual double unscaleX(double aX, double aZ) = 0;

	//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
	virtual double unscaleY(double aY, double aE/*eV*/, double aZ) = 0;

private:
	void InitRanges(double aZmax);
	CAutoDeletePtrArray<CTableFunction>	fData;
	CVector								fZ;
	double								fMaxDataZ;
	double								fMinE;
	double								fMaxE;
	double								fRangesZmax;
};

class CTableWithHeaderBackgr : public CTableWithHeaderReader, public IBackgroundSpectrum
{
public:
	CTableWithHeaderBackgr(const char* theFileName, char theCommentSimbol = '#'):
	  CTableWithHeaderReader(theFileName, theCommentSimbol){};

	/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z)=0;
	virtual bool init(){return read();};
};

/// This class was originally designed to read slightly modified Elmag background tables
/// (see script ../MCPropag/ElmagTest/Tables/makeTablesForMCPropag.sh)
/// The table file is supposed to use the following format
/// 	 z1				z2				z3				...			zMax
/// E1	 dN/dE(z1,E1)	dN/dE(z2,E1)	dN/dE(z3,E1)	...			dN/dE(zMax,E1)
/// E2	 dN/dE(z1,E2)	dN/dE(z2,E2)	dN/dE(z3,E2)	...			dN/dE(zMax,E2)
/// ................
/// Emax dN/dE(z1,Emax)	dN/dE(z2,Emax)	dN/dE(z3,Emax)	...			dN/dE(zMax,Emax)
///
/// where [E]=eV, [dN/dE] = 1/cm^3/eV
/// Linear interpolation in z, log(E) and log (N) is performed as in Elmag
/// By default for all z<z_min dN/dE(z) == dN/dE(zMin)
class MatrixBackground : public IBackgroundSpectrum{
public:
	MatrixBackground(string aTableFile, bool aIsComoving, bool aExtendToZero);

	virtual ~MatrixBackground();

	/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z);

	virtual double MaxZ() const;

	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;

	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;

private:
	CAutoDeletePtrArray<CVector>	fLogEdN_dE;
	CVector							fZ;
	CVector							fZindex;
	CDefaultTableFunc*				fZtoIndex;
	CVector							fLogE;
	CVector							fLogEindex;
	CDefaultTableFunc*				fLogEtoIndex;
	bool							fExtendToZero;
	bool							fIsComoving;
	char*							fBuffer;
	double							fMaxE;
	double							fMinE;
};

//plain table background (just file with two columns E/eV & E*dN/dE/cm^{-3}
//the background evolution can be implemented via optional aEvolution parameter
//which defines background with proportional evolution in this case the table defines background at z=0
class PlainTableBackground : public IBackgroundSpectrum
{
public:
	virtual double F(double E, double z);

	virtual double MaxZ() const;

	virtual double MaxE(double aZmax) const;

	virtual double MinE(double aZmax) const;

	virtual bool init();

	PlainTableBackground(std::string aTableFile, bool aConstComovingDensity, IBackgroundSpectrum* aEvolution=0, double aZmax=1e100);
private:
	bool fConstComovingDensity;
	double fZmax;
	SafePtr<CTableFunction>     fTable;
    SafePtr<IBackgroundSpectrum>  fEvolution;
};

#endif //#if !defined(TableBackground_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
