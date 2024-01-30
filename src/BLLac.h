/*
 * BLLac.h
 *
 *  Created on: Jun 5, 2015
 *      Author: ok
 */

#ifndef BLLAC_H_
#define BLLAC_H_
#include "TableFunc.h"

//BLLac luminosity-dependent density evolution (LDDE)
//parameterization form of Ueda, Y., Akiyama, M., Ohta, K., & Miyaji, T. 2003, ApJ, 598, 886
//with best fit parameters from arXiv/1311.5708v3
class BLLacLDDE : public IFunction{
public:
	BLLacLDDE();
	static double e(double z, double L/*in units of 10^48 erg/s*/);
	static double rho(double L/*in units of 10^48 erg/s*/);
	static void UnitTest();
	static double N(double z);//BLLac density (here used for self-test purposes only)
	static double NL(double z);// BLLac emission density
	//virtual double f(double z) const { return fTable->f(z); }//todo: tabulate function
	virtual double f(double z) const { return NL(z)/fNL0; }//normalize to 1 at z=0
private:
	//SafePtr<CLinearFunc> fTable;//todo: tabulate function
	static const double Lmin;
	static const double Lmax;
	double fNL0;//NL(0)
};

#endif /* BLLAC_H_ */
