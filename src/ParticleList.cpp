// ParticleFactory.cpp: implementation of the CParticleFactory class.
//
//////////////////////////////////////////////////////////////////////

#include "ParticleList.h"
#include "Parameters.h"
#include "Nucleus.h"

const int CParticleList::NucleiCount = EEndNuclei - EStartNuclei;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



CParticleList::CParticleList():
fMinNucleiA(2),//calculate nuclei propagation up to iMinNucleiA
fMaxNucleiA(NucleiCount+1),//calculate nuclei propagation starting from iMaxNucleiA
fNucleiEnabled(true)
{
	for(int particle=0;particle<(int)EEndAllParticles;particle++){
		fDoCalcSpectrum[particle] = particle<(int)EEndParticle;
		fExternal[particle] = false;
	}

	ReadSwitches();

	InitPropagatingParticlesIndex();
}

CParticleList::~CParticleList()
{

}

CParticleList* CParticleList::fInstance = 0;

void CParticleList::ReadSwitches()
{
	int particle=0;
	for(particle=0;particle<(int)EEndParticle;particle++)
	{
		if (particle>=EStartNuclei&&particle<EEndNuclei) {
			continue;//skip nuclei
		}
		const char* particleName = ParticleData::getParticleFileName((TParticle)particle);
		char key[64] = "calc_";
		strcat(key,particleName);
		Reader()->readBoolPar(key,fDoCalcSpectrum[particle]);
	}
	Reader()->readIntPar("MinNucleiA",fMinNucleiA);
	Reader()->readIntPar("MaxNucleiA",fMaxNucleiA);

	if (fMinNucleiA<2) fMinNucleiA=2;
	if (fMaxNucleiA>NucleiCount+1) fMaxNucleiA=NucleiCount+1;
	if (fMinNucleiA>fMaxNucleiA) {
		fMinNucleiA = fMaxNucleiA = 1;
		fNucleiEnabled = false;
	}

	for(particle=EStartNuclei;particle < EStartNuclei+fMinNucleiA-2;particle++)
		fDoCalcSpectrum[particle] = false;

	for(particle=EStartNuclei+fMaxNucleiA-1 ;particle < EEndNuclei;particle++)
		fDoCalcSpectrum[particle] = false;
}

void CParticleList::InitPropagatingParticlesIndex()
{
	if(!fPropagatingParticleIndex.length())
	{
		fPropagatingParticleIndex.create(EEndAllParticles);
		int i=0;
		for(int particle = 0; particle < EEndAllParticles; particle++)
		{
			if(IsEnabled((TParticle)particle))
			{
				fPropagatingParticleIndex[particle] = i++;
				fPropagatingParticles.add((TParticle)particle);
			}
			else
				fPropagatingParticleIndex[particle] = -1;
		}
	}
}


CVarArray<TParticle>*	ParticleIterator::fParticles = NULL;

ParticleIterator::ParticleIterator(TParticleIteratorType aType):
fPos(0),
fCurrentSubset(NULL)
{
	ASSERT(aType!=ENumberOfParticleIteratorTypes);
	if (fParticles==NULL) {
		init();
	}
	fCurrentSubset = fParticles + aType;
}

void ParticleIterator::init()
{
	fParticles = new CVarArray<TParticle>[ENumberOfParticleIteratorTypes];
	CParticleList* factory = CParticleList::Instance();
	int i;

	for(i=EStartLightParticle; i<EEndAllParticles; i++){
		if (factory->IsEnabled((TParticle)i)) {
			fParticles[EAll].add((TParticle)i);
		};
	}
	
	for(i=EStartLightParticle; i<EEndParticle; i++){
		if (factory->IsEnabled((TParticle)i)) {
			fParticles[EReal].add((TParticle)i);
		};
	}

	for(i=EStartEM;i<(int)EEndEM;i++){
		if (factory->IsEnabled((TParticle)i)) {
			fParticles[EEM].add((TParticle)i);
		};
	}
	if (factory->IsEnabled(EReserved5)) {
		fParticles[EEM].add(EReserved5);
	}
	
	for(i=(int)EEndEM;i<(int)EEndParticle;i++){
		if (factory->IsEnabled((TParticle)i)) {
			fParticles[ENonEM].add((TParticle)i);
		};
	}

	for(i=(int)EStartNuclei;i<(int)EEndNuclei;i++){
		if (factory->IsEnabled((TParticle)i)) {
			fParticles[ENuclei].add((TParticle)i);
		};
	}
}

