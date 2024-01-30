/*
 * Franceschini08EBL.cpp
 *
 *  Created on: Dec 22, 2013
 *      Author: ok
 */

#include "Franceschini08EBL.h"

static const char* files[] = {"0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0", 0};

Franceschini08EBL::Franceschini08EBL():
TableBackground("Franceschini08EBL", files, true)
{
}

double Franceschini08EBL::scaleX(double aE/*eV*/, double aZ)
{
	return log10(aE);
}

double Franceschini08EBL::unscaleX(double aX, double aZ)
{
	return pow(10., aX);
}

//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
double Franceschini08EBL::unscaleY(double aY, double aE/*eV*/, double aZ)
{
	double dl = (1.+aZ);
	return pow(10., aY)/(dl*dl*dl);//data is given in log scale in physical volume => converting to comoving
}
