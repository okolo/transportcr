//
// Custom analytical background
//

#ifndef PROPAGATION_CUSTOMBACKGROUND_H
#define PROPAGATION_CUSTOMBACKGROUND_H


#include "Background.h"

class CustomBackground :
        public IBackgroundSpectrum {
public:
    CustomBackground();
    virtual bool init();
    /* Photon spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
      (must be multiplied by (1+z)^3 before substituting farther)*/
    virtual double F(double E, double z);

    //returns maximal background red shift
    virtual double MaxZ() const;

    //returns maximal background energy in eV
    virtual double MaxE(double aZmax) const;

    //returns minimal background energy in eV
    virtual double MinE(double aZmax) const;
protected:
    double fEmin;
    double fEcut;
    double fNormN;
};


#endif //PROPAGATION_CUSTOMBACKGROUND_H
