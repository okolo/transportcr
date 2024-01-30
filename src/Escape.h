/*
 * Decay.h
 *
 *  Created on: May 13, 2013
 *      Author: ok
 */

#ifndef DECAY_H_
#define DECAY_H_

#include "Coupling.h"

namespace couplings {

class Escape: public Coupling {
	class EscapeChannel : public CouplingChannelT<Escape>
		{
		public:
			EscapeChannel(Escape* aCoupling, TParticle aPrim) : CouplingChannelT<Escape>(aCoupling, aPrim, aPrim){};
			virtual void Coef(CMatrixAddOnlyView& aCoef) const;
		};
public:
	Escape(TParticle aPrim);
private:

	//escape probability per unit time in internal units [T^-1]
	//for simple decay it would be R = 1/tau/(aE/m) , where tau is rest frame life time, m is particle mass
	virtual double EscapeRate(double aE) = 0;
};

class ConstEscape : public Escape
{
public:
	ConstEscape(TParticle aPrim, double aEscapeLengthMpc);
private:
	virtual double EscapeRate(double aE);
	double fRate;
};

/*
 * Escape of charged particles from volume filled with random magnetic field
 * **/
class MagneticDiffusionEscape : public Escape
{
public:
	MagneticDiffusionEscape(TParticle aPrim, double aEscapeLengthMpc, double aCriticalEnergy_eV);
private:
	virtual double EscapeRate(double aE);
	double fRate;
	double fCriticalEnergy;
};

/**
 * Escape of charged particles from a cone (limited by fixed angle) due to deflection by regular magnetic field
 * */
class MagneticDeflectionEscape : public Escape
{
public:
	MagneticDeflectionEscape(TParticle aPrim, double aBgauss, double aMaxDeflectionAngleDegrees);
private:
	virtual double EscapeRate(double aE);
	double fRateAtUnitE;
};

/**
 * Escape of charged particles from a beam directed from a source at z=Ranges().Zmax() to observer at z=0
 * due to deflection by regular magnetic field
 * */
class MagneticDeflectionJetEscape : public Escape
{
public:
	//at the time of introduction of this coefficient aEnergyScaleCorrection = 1.5
		//was needed to achieve good agreement with Monte Carlo simulations with l_cor=1Mpc, E=1e14 eV and B=1e-15 G
	MagneticDeflectionJetEscape(TParticle aPrim, double aBgauss, double aPSFhalfAngleDeg, double aJetOpenningAngleDeg, double aEnergyScaleCorrection=1.);
private:
	virtual double EscapeRate(double aE);
	double fRateAtUnitE;
	double fSourceD;
	double fsinPSF;
	double fcosPSF;
	double fsinJet;
};

} /* namespace couplings */
#endif /* DECAY_H_ */
