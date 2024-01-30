/*
 * Franceschini08EBL.h
 *
 *  Created on: Dec 22, 2013
 *      Author: ok
 */

#ifndef FRANCESCHINI08EBL_H_
#define FRANCESCHINI08EBL_H_

#include "TableBackground.h"

/// EBL from Table1,2 of http://arxiv.org/pdf/0805.1841v2.pdf
class Franceschini08EBL: public TableBackground {
public:
	Franceschini08EBL();
protected:
	//convert aE to internal x scale
	virtual double scaleX(double aE/*eV*/, double aZ);

	double unscaleX(double aX, double aZ);

	//convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
	virtual double unscaleY(double aY, double aE/*eV*/, double aZ);
};

#endif /* FRANCESCHINI08EBL_H_ */
