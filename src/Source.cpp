#include "Source.h"
#include "ParticleList.h"
#include "TimeZ.h"
#include "Parameters.h"
#include "Units.h"
#include "InjectionSpectra.h"

SourceTable::SourceTable(CInjectionSpectra* aInjection):
fWeight(1.),
fInjection(aInjection),
fT(-1),
fEnergyCorrectionFactor(1),
fEnergyDepWeight(0)
{
	create(EEndAllParticles,Ranges().nE());
}

double SourceTable::GetSource(TParticle _particle,
		double _E/* E/A for nuclei  */,
		int _i,
		double aEnergyCorrectionFactor,
		double _z) const
{
	if(fWeight==0)
		return 0.;
	double result = fInjection->ModifiedSpectrum(_particle, _i, _E*aEnergyCorrectionFactor*units.Eunit, _z);
	ASSERT_VALID_NO(result);

	return fWeight*BC().sF*_E*ParticleData::GetEnergyScaleFactor(_particle)*units.SpecUnit*result;
}

void SourceTable::calculateQ(double _t, double aEnergyCorrectionFactor)
{
	fT = _t;
	fEnergyCorrectionFactor = aEnergyCorrectionFactor;
	const CBinning& E = Ranges().midE();
	const int nn = Ranges().nE();
	double z = (_t>0.0) ? CTimeZ::t2z(_t) : 0.;
	double baseWeight = fWeight;
	for(int i=0; i<nn; i++) {
		double curE = E[i];
		if(fEnergyDepWeight)
			fWeight = baseWeight*fEnergyDepWeight->f(curE);
		FOR_ALL_REAL_PARTICLES_INVOLVED(particle) {
			(*this)[particle][i] = GetSource((TParticle) particle, curE, i, aEnergyCorrectionFactor, z);
		}
	}
	fWeight = baseWeight;
}

void SourceTable::SetWeight(double aWeight, IFunction* aWeightE)
{
	ASSERT_VALID_NO(aWeight);
	if(fEnergyDepWeight || aWeightE)
	{//always recalculate
		fWeight = aWeight;
		fEnergyDepWeight = aWeightE;
		if(aWeight==0.)
			reset();
		else
			calculateQ(fT, fEnergyCorrectionFactor);
	}
	else {
		if (fWeight > 0) {
			ASSERT_VALID_NO(aWeight / fWeight);
			multiply(aWeight / fWeight);
			fWeight = aWeight;
		}
		else if (aWeight > 0) {
			fWeight = aWeight;
			calculateQ(fT, fEnergyCorrectionFactor);
		}
	}
}

const double SourceTable::DefaultMomentSourceDt = 1.;//internal units

void SourceTable::MomentSource(Concentrations &aN, double aZ, bool aExternalOnly, IFunction* aEnergyDepWeight, double _dt) const
{
	const int nn = Ranges().nE();
	CParticleList* pl = CParticleList::Instance();

	double mem_weight=fWeight;
	IFunction* energyDepWeight = aEnergyDepWeight ? aEnergyDepWeight : fEnergyDepWeight;

	fWeight=_dt;
	const CBinning& E = Ranges().midE();

	for(int i=0;i<nn;i++)
	{
		double curE = E[i];
		if(energyDepWeight)
			fWeight = _dt*energyDepWeight->f(curE);
		FOR_ALL_REAL_PARTICLES_INVOLVED(particle){
			if(aExternalOnly && !pl->IsExternal(particle))
				continue;
			aN.SetBin((TParticle)particle, i, GetSource((TParticle)particle,curE,i,1.,aZ));
		}
	}

	fWeight=mem_weight;
}

