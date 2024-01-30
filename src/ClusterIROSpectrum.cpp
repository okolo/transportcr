// ClusterIROSpectrum.cpp: implementation of the CClusterIROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "ClusterIROSpectrum.h"
#include "Units.h"
#include <math.h>
#include "TableFunc.h"
#include "TableReader.h"
#include <fstream>

#define CLUSTER_DATA_DIR DATA_DIR "cluster" DIR_DELIMITER_STR

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const double CClusterIROSpectrum::Lsolar = 3.827e33;// erg/sec  (1 erg/s = 1e-7 W )
const double CClusterIROSpectrum::clusterSizeMpc = 3.0;//Mpc
const double CClusterIROSpectrum::numberOfGalaxies=1e3;

CClusterIROSpectrum::CClusterIROSpectrum():
iEllipticalGalaxiesFraction(NULL),
iEllipticalGalaxySpectrum(NULL),
iSpiralGalaxySpectrum(NULL),
iGalaxyConcentration(NULL),
iL(clusterSizeMpc*units.Mpc_cm /units.Lunit)
{
	
}

void CClusterIROSpectrum::printGalaxySpectra(ostream& aOut) const
{
	aOut << "#Total luminosity of elliptical(el) and spiral(sp) galaxy\n#E, eV\t\tdN_el(E)/dE x E\t\tdN_sp(E)/dE x E\n";
	//double sF = BC().sF;
	for(double E = fMinK_eV; E<=fMaxK_eV; E*=BC().ss1){
		double innerE = E*1e-6/units.Eunit;
		aOut << E*units.Eunit*1e6  << "\t\t" << iEllipticalGalaxySpectrum->f(innerE)/units.Tunit
			<< "\t\t" << iSpiralGalaxySpectrum->f(innerE)/units.Tunit << '\n';
	}
	aOut.flush();
}

bool CClusterIROSpectrum::init()
{
	//double sF = BC().sF;
	double backgrEmin = 0;
	double backgrEmax = 0;
	{
		CFilePtr fractionFile(FopenL(CLUSTER_DATA_DIR "fraction-elliptical", "rt"));
		CTableReader* reader = new CTableReader(fractionFile,2);
		CVarVector& x = reader->getColumn(0);
		int iMax = x.length();
		for(int i=0; i<iMax; i++)
			x[i] *= units.Mpc_cm /units.Lunit;

		iEllipticalGalaxiesFraction = new CLinearFunc(reader, 1.,0.);
	}
	{	
		CFilePtr ellipticalFile(FopenL(CLUSTER_DATA_DIR "elliptical", "rt"));
		CTableReader* ellipticalData = new CTableReader(ellipticalFile,2);
		CVarVector& x = ellipticalData->getColumn(0);
		CVarVector& y = ellipticalData->getColumn(1);
		int iMax = x.length();
		x.invert();
		y.invert();
		for(int i=0; i<iMax; i++)
		{//x=wave length in mkm
			double wavelength = x[i]*1e-4/units.Lunit;//internal units
			x[i] = 2.*Pi/wavelength;//photon energy in internal units
			double L = y[i]*1e10*Lsolar/units.MeVinESU/units.Eunit*units.Tunit;
			y[i] = L/x[i];//*sF;
		}
		iEllipticalGalaxySpectrum = new CDefaultTableFunc(ellipticalData);
		backgrEmin = x[0];
		backgrEmax = x[iMax-1];
	}
	{	
		CFilePtr spiralFile(FopenL(CLUSTER_DATA_DIR "spiral", "rt"));
		CTableReader* spiralData = new CTableReader(spiralFile,2);
		CVarVector& x = spiralData->getColumn(0);
		CVarVector& y = spiralData->getColumn(1);
		int iMax = x.length();
		for(int i=0; i<iMax; i++)
		{//x=frequency in Hz
			double frequency = x[i]*units.Tunit;//internal units
			x[i] = 2.*Pi*frequency;//photon energy in internal units
			// y = dL/dE * dE (watt)
			double L = y[i]*1e7/units.MeVinESU/units.Eunit*units.Tunit;
			y[i] = L/x[i];//*sF;// converting L to N(E)*E / time unit
		}
		iSpiralGalaxySpectrum = new CDefaultTableFunc(spiralData);
		if (x[0]<backgrEmin) 
			backgrEmin = x[0];
		if (x[iMax-1]>backgrEmax) 
			backgrEmax = x[iMax-1];

	}
	fMinK_eV = backgrEmin*units.Eunit*1e6;
	fMaxK_eV = backgrEmax*units.Eunit*1e6;
	/*
	iBackgroundMin = Ranges().backgroundE().findX(backgrEmin);
	if (iBackgroundMin<0) {
		throw "invalid IR/O background range (minimal limit)";
	}
	iBackgroundMax = Ranges().midBackgroundE().findX(backgrEmax);
	if (iBackgroundMax<0) {
		throw "invalid IR/O background range (maximal limit)";
	}
	*/
	iGalaxyConcentration = new TGalaxyConcentration(iL, numberOfGalaxies);

	const int accuracy = 200;
	TRayIntegral I1(iGalaxyConcentration, iEllipticalGalaxiesFraction, iR, iL );
	TRayIntegral I2(iGalaxyConcentration, &gUnitFunc, iR, iL );
	fSum1 = FuncUtils::integrate(-1.,1., &I1, &gUnitFunc, accuracy);
	fSum2 = FuncUtils::integrate(-1.,1., &I2, &gUnitFunc, accuracy);

	return true;
}

double CClusterIROSpectrum::F(double aE_eV, double z)
{
	double E = aE_eV*1e-6/units.Eunit;
	double fEl = iEllipticalGalaxySpectrum->f(E);
	double fSp = iSpiralGalaxySpectrum->f(E);
	return 0.5*((fEl-fSp)*fSum1 + fSp*fSum2);
}

//returns maximal background red shift
double CClusterIROSpectrum::MaxZ() const
{
	return 100;//evolution is not supported
}

//returns maximal background energy in eV
double CClusterIROSpectrum::MaxE(double aZmax) const
{
	return fMaxK_eV;
}

//returns minimal background energy in eV
double CClusterIROSpectrum::MinE(double aZmax) const
{
	return fMinK_eV;
}


CClusterIROSpectrum::~CClusterIROSpectrum()
{
	delete iEllipticalGalaxiesFraction;
	delete iEllipticalGalaxySpectrum;
	delete iGalaxyConcentration;
	delete iSpiralGalaxySpectrum;
}

void CClusterIROSpectrum::makeTests()
{
#define CHECK_SET_SIZE 6
	double checkRset[CHECK_SET_SIZE]={0, 0.1, 0.3, 0.5, 1, 2};
	try{
		CClusterIROSpectrum self;
		self.init();
		{
			string outFileName = plt_local_dir + "clusterIRO";
			ofstream out(outFileName.c_str());
			for(int i=0; i<CHECK_SET_SIZE; i++){
				self.setR(checkRset[i]*units.Mpc_cm /units.Lunit);
				CBackgroundTable bt(&self);
				bt.update(0);
				out << bt;
				out << "\n\n";
			}
			out.close();
		}
		{
			string outFileName = plt_local_dir + "LgalaxyIRO";
			ofstream out(outFileName.c_str());
			self.printGalaxySpectra(out);
			out.close();
		}

	}catch (const char* errMsg) {
		cerr << errMsg;
	}
}

const double TGalaxyConcentration::r1Mpc=0.01;//Mpc
const double TGalaxyConcentration::r2Mpc=0.25;//Mpc
const double TGalaxyConcentration::r3Mpc=1.0;//Mpc
const double TGalaxyConcentration::alpha1=0.51;
const double TGalaxyConcentration::alpha2=0.72;
const double TGalaxyConcentration::alpha3=0.58;

static double square(double aX){
	return aX*aX;
}

TGalaxyConcentration::TGalaxyConcentration(double aMaxR, double Ntot):
iNorm(1.),
r1(r1Mpc*units.Mpc_cm /units.Lunit),
r2(r2Mpc*units.Mpc_cm /units.Lunit),
r3(r3Mpc*units.Mpc_cm /units.Lunit),
iMaxR(aMaxR)
{
	CFunc funcSquare(square);
	double sum = funcSquare.integrate(0, iMaxR, this, 1000);
	iNorm = 0.25*Ntot/Pi/sum;
}

double TGalaxyConcentration::f(double aR/*Mpc*/) const{
	return iNorm/(pow(1.+aR/r1,alpha1)*pow(1.+aR/r2,alpha2)*pow(1.+aR/r3,alpha3));
	//return iNorm;//test
}

