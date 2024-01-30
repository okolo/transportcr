

#include "ParticleData.h"
#include "Addfunc.h"
#include "FilePtr.h"
#include <math.h>
#include "main.h"
#include "TableFunc.h"
#include "TableReader.h"
#include "Concentrations.h"
#include "Nucleus.h"
#include "Units.h"

const double ParticleData::particleMassesMeV[]={
	0.51099906,//EElectron
	0.51099906,//EPositron
	0.,//EPhoton
	0.,//ENeutrinoE
	0.,//ENeutrinoM
	0.,//ENeutrinoT
	0.,//ENeutrinoAE
	0.,//ENeutrinoAM
	0.,//ENeutrinoAT
	939.56563,//ENeutron
	938,27231,//EProton
	-1.//end
};

const char*	ParticleData::defaultParticleFileNames[EEndAllParticles] = {
	//elementary particles
	"e",
	"ae",
	"ph",
	"en",
	"mn",
	"tn",
	"aen",
	"amn",
	"atn",
	"n",
	"p",
	//nuclei
	"A02",
	"A03",
	"A04",
	"A05",
	"A06",
	"A07",
	"A08",
	"A09",
	"A10",
	"A11",
	"A12",
	"A13",
	"A14",
	"A15",
	"A16",
	"A17",
	"A18",
	"A19",
	"A20",
	"A21",
	"A22",
	"A23",
	"A24",
	"A25",
	"A26",
	"A27",
	"A28",
	"A29",
	"A30",
	"A31",
	"A32",
	"A33",
	"A34",
	"A35",
	"A36",
	"A37",
	"A38",
	"A39",
	"A40",
	"A41",
	"A42",
	"A43",
	"A44",
	"A45",
	"A46",
	"A47",
	"A48",
	"A49",
	"A50",
	"A51",
	"A52",
	"A53",
	"A54",
	"A55",
	"A56",
	//reserved for future (to preserve output spectrum file format)
	"reserved1",
	"reserved2",
	"reserved3",
	"reserved4",
	"reserved5",
	"reserved6"
};

TParticle ParticleData::Antiparticle2Particle(const TParticle _antiparticle)
{
	TParticle result(_antiparticle);
	switch((int)_antiparticle)
	{
	case EPositron:
		result = EElectron;
		break;
	case ENeutrinoAE:
		result = ENeutrinoE;
		break;
	case ENeutrinoAM:
		result = ENeutrinoM;
		break;
	case ENeutrinoAT:
		result = ENeutrinoT;
		break;
	default:
		result = _antiparticle;
	}
	return result;
}

TParticle ParticleData::GetAntiparticle(const TParticle _particle)
{
	TParticle result(_particle);
	switch((int)_particle)
	{
	case EPositron:
		result = EElectron;
		break;
	case ENeutrinoAE:
		result = ENeutrinoE;
		break;
	case ENeutrinoAM:
		result = ENeutrinoM;
		break;
	case ENeutrinoAT:
		result = ENeutrinoT;
		break;
	case EElectron:
		result = EPositron;
		break;
	case ENeutrinoE:
		result = ENeutrinoAE;
		break;
	case ENeutrinoM:
		result = ENeutrinoAM;
		break;
	case ENeutrinoT:
		result = ENeutrinoAT;
		break;
	default:
		result = _particle;
	}
	return result;
}

const  double ParticleData::meanNucleonMassMeV=939.0; /* mean nuclon mass (MeV) */

const char* ParticleData::getParticleFileName(const TParticle aParticle)
{
	return defaultParticleFileNames[aParticle];
}

double ParticleData::getParticleMass(const TParticle _particle)
{
	if(_particle<EStartNuclei)
		return particleMassesMeV[_particle]/units.Eunit;
	if(_particle < EEndNuclei)
		return meanNucleonMassMeV*CNucleus::getA(_particle)/units.Eunit;
	ThrowError("CParticle::getParticleMass unsupported argument");
	return 0.;
}

int ParticleData::getElectricCharge(const TParticle _particle)
{
	switch((int)_particle)
	{
	case EPositron:
	case EProton:
		return 1;
	case ENeutrinoAE:
	case ENeutrinoE:
	case ENeutrinoAM:
	case ENeutrinoM:
	case ENeutrinoAT:
	case ENeutrinoT:
	case EPhoton:
	case ENeutron:
		return 0;
	case EElectron:
		return -1;
	}

	if(_particle < EEndNuclei)
		return CNucleus::getZ(_particle);
	ThrowError(ToString("getElectricCharge() is not implemented for argument ") + getParticleFileName(_particle));
	return 0;
}

double ParticleData::GetEnergyScaleFactor(TParticle aParticle)
{
	// Qi coef need to be adjusted for nuclei due to different energy binning, namely
	// since Qi = deltaEi * Q where delteEi = BC().sF * Enuclei = BC().sF * A * Ei
	return (aParticle<EEndNuclei && aParticle>=EStartNuclei) ? CNucleus::getA(aParticle) : 1;
}
