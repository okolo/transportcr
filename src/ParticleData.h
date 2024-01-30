#if !defined(PARTICLEDATA_H_INCLUDED)
#define PARTICLEDATA_H_INCLUDED

#include "Vector.h"
#include "TParticle.h"


class ParticleData
{
public:
	static double GetEnergyScaleFactor(TParticle aParticle);
	static TParticle Antiparticle2Particle(const TParticle _antiparticle);
	static TParticle GetAntiparticle(const TParticle _particle);
	static const char* getParticleFileName(const TParticle _particle);
	static int getElectricCharge(const TParticle _particle);
	static double getParticleMass(const TParticle _particle);
	static const  double meanNucleonMassMeV;
private:
	static const double particleMassesMeV[];
	static const char*	defaultParticleFileNames[EEndAllParticles];
};

#define FOR_ALL_PARTICLES(particleIndex) for(int particleIndex=0;particleIndex<(int)EEndAllParticles;particleIndex++)
#define FOR_ALL_REAL_PARTICLES(particleIndex) for(int particleIndex=0;particleIndex<(int)EEndParticle;particleIndex++)

#endif// end of file
