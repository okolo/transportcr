// ClusterIROSpectrum.h: interface for the CClusterIROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GALACTICIROSPECTRUM_H__5BB8B716_8E24_417C_95C1_CD6CB25C2396__INCLUDED_)
#define AFX_GALACTICIROSPECTRUM_H__5BB8B716_8E24_417C_95C1_CD6CB25C2396__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Function.h"
#include "Background.h"
#include <math.h>

class CClusterIROSpectrum : public IBackgroundSpectrum
{
public:
	CClusterIROSpectrum();
	virtual ~CClusterIROSpectrum();
	inline void setR(double aR){iR=aR;};//[aR]=internal units

	static void makeTests();
	void printGalaxySpectra(ostream& aOut) const;

	virtual bool init();

	virtual double F(double E, double z);

	//returns maximal background red shift
	virtual double MaxZ() const;

	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;

	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;


protected:
	void update(double /*aZ*/);
private:
	
	IFunction*  iEllipticalGalaxiesFraction;
	IFunction*  iEllipticalGalaxySpectrum;
	IFunction*  iSpiralGalaxySpectrum;
	IFunction*	iGalaxyConcentration;
	double      iNorm;
	double		iL;//cluster size, internal units
	double		iR;//current distance from the cluster center, internal units
	static const double Lsolar;//solar luminosity, erg/sec
	static const double clusterSizeMpc;//cluster size, Mpc
	static const double numberOfGalaxies;
	//int			iBackgroundMin;
	//int			iBackgroundMax;
	double		fMinK_eV;
	double		fMaxK_eV;
	double		fSum1;
	double		fSum2;

class TRayFunction : public IFunction{
	private:
		
		IFunction*  iConcentration;
		IFunction*  iWeight;
		double		iMu;
		double		iR;
		double		iL;
	public:
		TRayFunction(IFunction* aConcentration, //concentration
			IFunction* aWeight,	//Weight function
			double aR, //current distance from the center
			double aL  //max distance from the center (size of cluster)
			):iConcentration(aConcentration),iWeight(aWeight),iMu(-2),iR(aR),iL(aL){};

		double f(double r) const{
			double rho = sqrt(iR*iR+r*r-2.*iR*r*iMu);
			return iConcentration->f(rho) * iWeight->f(rho);
		};
		
		inline double setMu(double aMu)//returns rMax
		{
			iMu=aMu;
			return iR*iMu+sqrt(iR*iMu*iR*iMu+iL*iL-iR*iR);
		};
	};
	
	class TRayIntegral : public IFunction{
	public:
		TRayIntegral(IFunction* aConcentration, //concentration
			IFunction* aWeight, //Weight function
			double aR, //current distance from the center
			double aL):iBase(aConcentration, aWeight, aR, aL ),iAccuracyCoef(50./aL){};
		double f(double aMu/*cos(Theta)*/) const{
			double rMax = ((TRayFunction&)iBase).setMu(aMu);
			int ac = (int)(rMax*iAccuracyCoef);
			ac *= 2;
			return ac?FuncUtils::integrate(0,rMax,(TRayFunction*)(&iBase),&gUnitFunc,ac):0;
		}
	private:
		TRayFunction	 iBase;
		double			 iAccuracyCoef;
	};
};

class TGalaxyConcentration : public IFunction{
public:
	TGalaxyConcentration(double aMaxR, double Ntot=1e3);
	virtual double f(double _x) const;
	inline double maxR(){return iMaxR;};
private:
	static const double r1Mpc;
	static const double r2Mpc;
	static const double r3Mpc;
	double iNorm;
	double r1;
	double r2;
	double r3;

	static const double alpha1;
	static const double alpha2;
	static const double alpha3;
	double iMaxR;
};

#endif // !defined(AFX_GALACTICIROSPECTRUM_H__5BB8B716_8E24_417C_95C1_CD6CB25C2396__INCLUDED_)
