/*
 * Galaxy.h
 *
 *  Created on: Aug 4, 2014
 *      Author: mac
 */

#ifndef GALAXY_H_
#define GALAXY_H_

#include "PropagEngine.h"
#include "InjectionSpectra.h"
#include "Medium.h"

///Rectangular coordinate system (CS) with center in Milky Way center, (x,y) plane parallel to disk
///and Sun's position given by (-R_sun,0,-Z_sun)
class DiskCenterCS {
public:
	DiskCenterCS();

	//Convert from galactic coordinate system to disk center CS
	//l,b should be given in degrees
	//output values (x,y,z) are given in kpc
	void GalacticToDiskCenterCS(double l, double b, double rKpc, double &x, double &y, double &z) const;

	//Calculate distance from galactic center using galactic coordinates l,b and distance from sun rKpc in kpc
	double GetDistanceFromCenter(double l, double b, double rKpc) const;

	//Calculate distance from sun to object in the direction l,b in galactic coordinates at given distance from galactic center
	//l,b should be given in degrees
	//return value is given in kpc
	//if aDistanceFromDiskCenterKpc<R_sun two solutions exist
	//the largest distance is returned by this function
	double GetDistanceFromSun(double l, double b, double aDistanceFromDiskCenterKpc) const;
	static const double Z_sun;
	static const double R_sun;
protected:
	const double sinBeta;//CO rotation angle sin for convertion to/from galactic coordinates
	const double cosBeta;//CO rotation angle cos for convertion to/from galactic coordinates
};

class GalacticDensityDistribution
{
public:
	// return true if this distribution is isotropic in galactic center coordinate system
	// i.e. density depends only on distance from galactic center
	virtual bool IsIsotropic() const = 0;
	virtual double MaxDistanceKpc(double l, double b) const = 0;
	virtual double Density(double l, double b, double rKpc) const = 0;
	virtual IFunction* EnergyWeight(double l, double b, double rKpc) const { return 0; }
	virtual ~GalacticDensityDistribution();
	//Get external spectra if it is used
	virtual void GetExternalSpectra(Concentrations& aN1, double l, double b, double rKpc, const SourceTable* aSource) const {}
	static GalacticDensityDistribution* Create(IParReader* aReader, double aSmoothSizeKpc);

	//Return true if jacobian recalculation is needed
	virtual bool UpdateMedium(Medium& aMedium, double l, double b, double rKpc) const { return false; }
	virtual Medium* CreateMedium() const { return new Medium(); }
};

/// Calculate directional or average flux from galactic source distribution at Sun's position
/// The source distribution is given as function of of galactic coordinates (l,b)
/// and distance from the sun r
class GalaxyPropagEngine : protected CPropagEngine, CouplingParameters
{
public:
	GalaxyPropagEngine(CInjectionSpectra* aInjection, double aMicroStep);
	~GalaxyPropagEngine();
	void run(Concentrations& aN1, double l, double b);//calculate flux from direction l, b in Galactic coordinates (in degrees)
	//calculate total flux (l,b are given in degrees)
	void run(Concentrations& aN1, double lMin, double lMax, double bMin, double bMax);

protected:
	virtual bool beforeStep(int stepNo, double t, double dt);
	virtual bool start();
private:
	TPropagCoef									fPropagCoef;
	SmartPtr<CInjectionSpectra> 				fInjection;
	SafePtr<GalacticDensityDistribution>	fDensity;
	double										fCurL;
	double										fCurB;
	double										fAngleStep;
	SafePtr<IFunction>							fCurWeightE;
    bool fVerboseOutput;
};

//calculate total flux of halo as an external object in [Tunit]^-1
//current implementation is oversimplified
//it is only valid for ISOTROPIC source of NONINTERACTING particles
//TODO: build proper implementation
class GalaxyEffectiveSourceEngine : protected CPropagEngine, CouplingParameters
{
public:
	GalaxyEffectiveSourceEngine(CInjectionSpectra* aInjection, double aMicroStep);
	~GalaxyEffectiveSourceEngine();

	//calculate total flux of halo as an external object in [Tunit]^-1
	void run(Concentrations& aN1);

protected:
	virtual bool beforeStep(int stepNo, double t, double dt);
	virtual bool start();
private:
	TPropagCoef									fPropagCoef;
	SmartPtr<CInjectionSpectra> 				fInjection;
	SafePtr<GalacticDensityDistribution>		fDensity;
	SafePtr<IFunction>							fCurWeightE;
	double										fR;
};

#endif /* GALAXY_H_ */
