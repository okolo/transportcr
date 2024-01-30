/*
 * CustomCouplings.cpp
 *
 *  Created on: May 14, 2013
 *      Author: ok
 */

#include "Coupling.h"
#include "Escape.h"
#include "Parameters.h"

using namespace couplings;

class CouplingListCustomPars : Parameters
{
public:
	static double Bgauss;
	static double MaxDeflectionDeg;
	static double JetOpenningAngleDeg;
	static double DeflectionEnergyScaleCorrection;
	CouplingListCustomPars()
	{
		if(ready) return;
		Reader()->readDoublePar("ccBgauss", Bgauss);
		Reader()->readDoublePar("ccMaxDeflectionDegrees", MaxDeflectionDeg);
		Reader()->readDoublePar("ccJetOpenningAngleDeg", JetOpenningAngleDeg);
		Reader()->readDoublePar("ccDeflectionEnergyScaleCorrection", DeflectionEnergyScaleCorrection);
		ready = true;
	}
private:
	static bool ready;
};

double CouplingListCustomPars::Bgauss = 1e-15;
double CouplingListCustomPars::MaxDeflectionDeg = 0.03;
double CouplingListCustomPars::JetOpenningAngleDeg = 180.;
double CouplingListCustomPars::DeflectionEnergyScaleCorrection = 1.5;
bool CouplingListCustomPars::ready = false;

void CouplingList::CustomInit()
{
	CouplingListCustomPars customCouplingParams;

	double Bgauss=customCouplingParams.Bgauss;
	double angleDegrees=customCouplingParams.MaxDeflectionDeg;
	double angleJet=customCouplingParams.JetOpenningAngleDeg;
	//FOR_ALL_PARTICLES_INVOLVED(particle)
	//{
	//	if(getElectricCharge(particle))
	//		fCouplings.add(new MagneticDeflectionEscape(particle, Bgauss, angleDegrees));
	//}
	//only electrons and positrons are deflected by correlated B
	//for protons correlation length is smaller than energy loss length and therefore diffusion mechanism should be used
	if(Bgauss<=0.)
		return;
	if(angleJet<=0.)
	{
		fCouplings.add(new MagneticDeflectionEscape(EElectron, Bgauss, angleDegrees));
		fCouplings.add(new MagneticDeflectionEscape(EPositron, Bgauss, angleDegrees));
	}
	else
	{
		fCouplings.add(new MagneticDeflectionJetEscape(EElectron, Bgauss, angleDegrees, angleJet,customCouplingParams.DeflectionEnergyScaleCorrection));
		fCouplings.add(new MagneticDeflectionJetEscape(EPositron, Bgauss, angleDegrees, angleJet,customCouplingParams.DeflectionEnergyScaleCorrection));
	}
}
