/*
 * Sarkar1005Background.h
 *
 *  Created on: May 19, 2011
 *      Author: ok
 */

#ifndef SARKAR1005BACKGROUND_H_
#define SARKAR1005BACKGROUND_H_

#include "Background.h"
#include "TableFunc.h"

enum SarkarBackgrType
{
	SarkarRadio=0,
	SarkarEBL=1
};

class Sarkar1005Background: public IBackgroundSpectrum
{
public:
	Sarkar1005Background(SarkarBackgrType aType);
	virtual ~Sarkar1005Background();
	virtual double F(double E, double z);
	virtual double MaxZ() const { return fMaxZ; }
	virtual double MaxE(double aZmax) const;
	virtual double MinE(double aZmax) const;
private:
	CTableFunction* f0;//spectrum at z=0
	IFunction* fEvol;//evolution func
	double	   fMaxZ;
};

#endif /* SARKAR1005BACKGROUND_H_ */
