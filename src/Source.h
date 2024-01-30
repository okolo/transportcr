
#ifndef SOURCE_H_
#define SOURCE_H_

#include "Vector.h"
#include "Concentrations.h"
#include "InjectionSpectra.h"

class SourceTable : public CMatrix{
public:
	SourceTable(CInjectionSpectra* aInjection);

protected:
	double GetSource(TParticle _particle,
			double _E,// E/A for nuclei
			int _i,
			double aEnergyCorrectionFactor,
			double _z) const;
public:
	//void GetEnergyReceived(class Concentrations& aN,double _dt);
	void MomentSource(Concentrations& aN, double aZ, bool aExternalOnly, IFunction* aEnergyDepWeight = 0, double _dt = DefaultMomentSourceDt) const;
	void calculateQ(double _t=-1.0, double aEnergyCorrectionFactor = 1.);
	void SetWeight(double aWeight, IFunction* aWeightE=0);
	static const double DefaultMomentSourceDt;

private:
	mutable double fWeight;
	double fT;
	double fEnergyCorrectionFactor;
	SmartPtr<CInjectionSpectra> fInjection;
	IFunction* fEnergyDepWeight;
};

#endif /* SOURCE_H_ */
