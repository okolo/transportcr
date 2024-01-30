// ParticleFactory.h: interface for the CParticleFactory class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARTICLEFACTORY_H__1A33E7C6_5112_4D0E_9930_0E0CB893C25A__INCLUDED_)
#define AFX_PARTICLEFACTORY_H__1A33E7C6_5112_4D0E_9930_0E0CB893C25A__INCLUDED_


#include "Vector.h"
#include "TParticle.h"
#include "Parameters.h"

class IParReader;
class IParWriter;

class CParticleList : protected Parameters
{
public:
	CParticleList();
	virtual ~CParticleList();

	inline static CParticleList* Instance(){
		if(!fInstance)
			fInstance = new CParticleList();
		return fInstance;
	};

	inline static CParticleList* Set(CParticleList* aList){
		CParticleList* result = fInstance;
		fInstance = aList;
		return result;
	};

	inline const Array<int, int&>& GetPropagatingParticleIndex() const
	{
		return fPropagatingParticleIndex;
	}

	inline const Array<TParticle, TParticle&>& GetPropagatingParticles() const
	{
		return fPropagatingParticles;
	}

	inline int NumberOfPropagatingParticles() const
	{
		return fPropagatingParticles.length();
	}

	//if external flag is set, particle spectrum is not modified in Jacobian::FastEvolve and CPropagEngine::step
	inline bool 				IsExternal(TParticle aParticle) const {return fExternal[aParticle];}
	inline void					SetExternal(TParticle aParticle, bool aExternal=true){ fExternal[aParticle] = aExternal;}
	inline bool					IsEnabled(TParticle aParticle) const {return fDoCalcSpectrum[aParticle];}
	inline bool					NucleiEnabled() const {return fNucleiEnabled;}
	inline int					MaxNucleiA() const {return fMaxNucleiA;}
	static const int			NucleiCount;
private:
	void ReadSwitches();
	void InitPropagatingParticlesIndex();

	bool						fDoCalcSpectrum[EEndAllParticles];
	bool						fExternal[EEndAllParticles];
	int							fMinNucleiA;//calculate nuclei propagation up to iMinNucleiA
	int							fMaxNucleiA;//calculate nuclei propagation starting from iMaxNucleiA
	bool						fNucleiEnabled;//is set automatically depending on iMinNucleiA & iMaxNucleiA

	static CParticleList*		fInstance;
	//particle index used in matrix calculations (includes only propagating particles)
	Array<int, int&>			fPropagatingParticleIndex;
	CVarArray<TParticle>		fPropagatingParticles;
};

enum TParticleIteratorType{
		EAll = 0,
		EReal,
		EEM,
		ENonEM,
		ENuclei,
		ENumberOfParticleIteratorTypes
	};

class ParticleIterator{
public:
	
	ParticleIterator(TParticleIteratorType aType=EAll);
	inline int numberOfParticles() const {return fCurrentSubset->length();}
	inline TParticle nextParticle(){return fPos<fCurrentSubset->length()?(*fCurrentSubset)[fPos++]:EEndAllParticles;};
	inline bool hasMoreParticles(){return fPos<fCurrentSubset->length();};
	
private:
	static void	init();
	int								fPos;
	static CVarArray<TParticle>*	fParticles;
	CVarArray<TParticle>*			fCurrentSubset;
};

#define FOR_ALL_PARTICLES_INVOLVED_T(particle, type) ParticleIterator it_##particle(type);\
for(TParticle particle = it_##particle.nextParticle(); particle!=EEndAllParticles; particle = it_##particle.nextParticle())

#define FOR_ALL_PARTICLES_INVOLVED(particle) FOR_ALL_PARTICLES_INVOLVED_T(particle, EAll)
#define FOR_ALL_REAL_PARTICLES_INVOLVED(particle) FOR_ALL_PARTICLES_INVOLVED_T(particle, EReal)
#define FOR_EM_CASCADE_INVOLVED(particle) FOR_ALL_PARTICLES_INVOLVED_T(particle, EEM)
#define FOR_NON_EM_CASCADE_INVOLVED(particle) FOR_ALL_PARTICLES_INVOLVED_T(particle, ENonEM)
#define FOR_ALL_NUCLEI_INVOLVED(particle) FOR_ALL_PARTICLES_INVOLVED_T(particle, ENuclei)


#endif // !defined(AFX_PARTICLEFACTORY_H__1A33E7C6_5112_4D0E_9930_0E0CB893C25A__INCLUDED_)
