/*
 * Decay.cpp
 *
 *  Created on: May 13, 2013
 *      Author: ok
 */

#include <math.h>
#include "Escape.h"
#include "Units.h"
#include "TimeZ.h"
#include "Ranges.h"
#include "ParticleData.h"

namespace couplings {

Escape::Escape(TParticle aPrim)
{
	AddChannel(new EscapeChannel(this, aPrim));
}

void Escape::EscapeChannel::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CBinning& energy = Ranges().midE();

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, -fCoupling.EscapeRate(energy[iPrim]));
	}
}

ConstEscape::ConstEscape(TParticle aPrim, double aEscapeLengthMpc):
		Escape(aPrim)
{
	fRate = 1./(aEscapeLengthMpc*units.Mpc_cm /units.Lunit);
}

double ConstEscape::EscapeRate(double aE)
{
	return fRate;
}

MagneticDiffusionEscape::MagneticDiffusionEscape(TParticle aPrim, double aEscapeLengthMpc, double aCriticalEnergy_eV):
		Escape(aPrim)
{
	fRate = 1./(aEscapeLengthMpc*units.Mpc_cm /units.Lunit);
	fCriticalEnergy = aCriticalEnergy_eV*1e-6/units.Eunit;
}

double MagneticDiffusionEscape::EscapeRate(double aE)
{
	return aE>fCriticalEnergy ? fRate : fRate*aE/fCriticalEnergy;
}

MagneticDeflectionEscape::MagneticDeflectionEscape(TParticle aPrim, double aBgauss, double aMaxDeflectionAngleDegrees):
		Escape(aPrim)
{
	int q = abs(ParticleData::getElectricCharge(aPrim));
	if(q==0)
	{
		fRateAtUnitE = 0.;
		return;
	}
	fRateAtUnitE = 5.9157134470e-15*aBgauss/units.Eunit/(aMaxDeflectionAngleDegrees/180.*M_PI)/units.Eunit*q/exp(1.);
}

double MagneticDeflectionEscape::EscapeRate(double aE)
{
	return fRateAtUnitE/aE;
}

MagneticDeflectionJetEscape::MagneticDeflectionJetEscape(TParticle aPrim, double aBgauss, double aPSFhalfAngleDeg, double aJetOpenningAngleDeg, double aEnergyScaleCorrection):
		Escape(aPrim),
		fRateAtUnitE(0),
		fSourceD(-1.),//will be calculated when needed
		fsinPSF(sin(aPSFhalfAngleDeg*M_PI/180.)),
		fcosPSF(cos(aPSFhalfAngleDeg*M_PI/180.)),
		fsinJet( aJetOpenningAngleDeg<90 ? sin(aJetOpenningAngleDeg*M_PI/180.) : 1.)
{
	int q = abs(ParticleData::getElectricCharge(aPrim));
	if(q==0 || aPSFhalfAngleDeg>=90.)
		return;
	fRateAtUnitE = 5.9157134470e-15*aBgauss/units.Eunit/units.Eunit*q/exp(1.)/aEnergyScaleCorrection;
}

double MagneticDeflectionJetEscape::EscapeRate(double aE)
{
	if(fSourceD<0.)
		fSourceD = CTimeZ::z2d(Ranges().Zmax());
	double d = CTimeZ::z2d(redshift.z());
	double r = sqrt(fSourceD*fSourceD+d*d-2.*fSourceD*d*fcosPSF);//distance to particle at given z provided that it arrives to observer at the edge of PSF
	double maxDeflectionSinByPSF = fSourceD/r*fsinPSF;
	double maxDeflectionSinByJet = fSourceD/d*fsinJet;
	double maxDeflectionSin = maxDeflectionSinByPSF>maxDeflectionSinByJet ? maxDeflectionSinByJet : maxDeflectionSinByPSF;//pick stronger condition
	double maxDeflectionAngle = asin(maxDeflectionSin);// maximal deflection angle (in radians) allowed at current distance from the observer
	return fRateAtUnitE/aE/maxDeflectionAngle;
}

} /* namespace couplings */
