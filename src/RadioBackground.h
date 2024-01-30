
#ifndef RADIOBACKGROUND_H_
#define RADIOBACKGROUND_H_

#include "Background.h"

//Astronomy & Astrophysics 365, 294-300 (2001) Calibration of low-frequency radio telescopes
//using the galactic background radiation (G.A. Dulk , ...)
//currently evolution with z is not supported
class DulkRadioBackground: public IBackgroundSpectrum {
public:
	DulkRadioBackground(bool aIncludeGalactic, bool aIsComovingConst = false);

	/* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z);

	//returns maximal background red shift
	virtual double MaxZ() const;

	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;

	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;
private:
	const bool 	includeGalactic;
	const double minFreq;
	const double maxFreq;
	static const double eVtoMHz;
	bool IsComovingConst;
};

//radio spectrum at z=0 : R.J. Protheroe, P.L. Biermann  astro-ph/9605119 (solid curves from Figure 5)
class ProtheroeRadioBackground: public IBackgroundSpectrum {
public:
	ProtheroeRadioBackground(bool aUpper, bool aIsComovingConst = false);

	/* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
	virtual double F(double E, double z);

	//returns maximal background red shift
	virtual double MaxZ() const;

	//returns maximal background energy in eV
	virtual double MaxE(double aZmax) const;

	//returns minimal background energy in eV
	virtual double MinE(double aZmax) const;
private:
	static const int dataSize;
	static const double x[];
	static const double y_low[];
	static const double y_up[];
	const double* y;
	bool IsComovingConst;
};


#endif /* RADIOBACKGROUND_H_ */
