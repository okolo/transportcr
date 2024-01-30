//
// Created by ok on 01/09/2016.
//

#ifndef PROPAGATION_STECKER16BACKGROUND_H
#define PROPAGATION_STECKER16BACKGROUND_H

#include "Background.h"
#include "Vector.h"
#include "TableFunc.h"
#include <string>
#include <vector>
#include "DataReader.h"
using namespace std;

class Stecker16Background : public IBackgroundSpectrum{
public:
    Stecker16Background(string aTableFile):fMatrixFunction(aTableFile){};

    virtual ~Stecker16Background(){};

    /* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
      (must be multiplied by (1+z)^3 before substituting farther)*/
    virtual double F(double E, double z);

    virtual double MaxZ() const { return fMatrixFunction.MaxArg(1);}

    //returns maximal background energy in eV
    virtual double MaxE(double aZmax) const { return fMatrixFunction.MaxArg(2);}

    //returns minimal background energy in eV
    virtual double MinE(double aZmax) const { return fMatrixFunction.MinArg(2);}

private:
    MatrixFunction  fMatrixFunction;
};

class Stecker16LowerBackground : public Stecker16Background{
public:
    Stecker16LowerBackground();
};

class Stecker16UpperBackground : public Stecker16Background{
public:
    Stecker16UpperBackground();
};
#endif //PROPAGATION_STECKER16BACKGROUND_H
