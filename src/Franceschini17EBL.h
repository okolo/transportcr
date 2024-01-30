//
// Created by ok on 29/05/18.
//

#ifndef PROPAGATION_FRANCESCHINI17EBL_H
#define PROPAGATION_FRANCESCHINI17EBL_H


#include "TableBackground.h"

class Franceschini17EBL :
        public TableBackground
{
public:
    Franceschini17EBL();
protected:
    //convert aE to internal x scale
    virtual double scaleX(double aE/*eV*/, double aZ);

    double unscaleX(double aX, double aZ);

    //convert internal y scale to output spectrum E*dn/dE in sm^-3 in comoving volume
    virtual double unscaleY(double aY, double aE/*eV*/, double aZ);
};


#endif //PROPAGATION_FRANCESCHINI17EBL_H
