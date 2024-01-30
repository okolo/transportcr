#if !defined(Kneiske1001IROSpectrum_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
#define Kneiske1001IROSpectrum_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_
#include "TableBackground.h"

class Kneiske1001IROSpectrum :
	public TableBackground
{
public:
	Kneiske1001IROSpectrum(void);
protected:
	//convert aE to internal x scale
	virtual double scaleX(double aE/*eV*/, double aZ);

	//convert internal x scale to energy in eV
	virtual double unscaleX(double aX, double aZ);

	//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
	virtual double unscaleY(double aY, double aE/*eV*/, double aZ);
};

class ElmagKneiskeMinimal : public MatrixBackground
{
public:
	/// Kneiske Lower limit Background used in Elmag
	/// There is an error in Elmag 1.02 and earlier, the tables are treated as physical concentrations
	/// But in fact they are built for comoving frame
	/// The parameter here is left just for comparison to Elmag 1.02
	/// In Elmag 2.0 this bug was fixed
	ElmagKneiskeMinimal(bool aIsComoving = true);
};

#endif //#if !defined(Kneiske1001IROSpectrum_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
