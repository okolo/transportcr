//
// Created by ok on 5/7/23.
//

#include "CustomInjSpectrum.h"
#include "TableFunc.h"

CCustomInjSpectrum::CCustomInjSpectrum()
{

}

void CCustomInjSpectrum::init()
{

}

double E[] = {0.5e6, 18e6, 251e6};
CVector vE(E, 3);
double E2F[] = {0.32, 0.9, 692};
CVector vE2F(E2F, 3);

PropagDepreciated::CLogScaleFunc specE2(vE, vE2F);
//CLinearFunc specE2(vE, vE2F);

double CCustomInjSpectrum::Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ)
{// https://www.overleaf.com/project/6453d74335a88241cb2f84b3
    return specE2(aE)/aE/aE;
}
