#if !defined(TPARTICLE_H_INCLUDED)
#define TPARTICLE_H_INCLUDED

enum _TParticle
{
	EStartLightParticle = 0,
	EStartEM  = EStartLightParticle,
	EElectron = EStartEM,
	EPositron,
	EPhoton,
	EStartNeutrino,
	ENeutrinoE = EStartNeutrino,
	ENeutrinoM,
	ENeutrinoT,
	ENeutrinoAE,
	ENeutrinoAM,
	ENeutrinoAT,
	EEndNeutrino,
	ENeutron = EEndNeutrino,
	EProton,
	EEndLightParticle,
	EStartNuclei = EEndLightParticle,
	EA02 = EStartNuclei,
	EA03,
	EA04,
	EA05,
	EA06,
	EA07,
	EA08,
	EA09,
	EA10,
	EA11,
	EA12,
	EA13,
	EA14,
	EA15,
	EA16,
	EA17,
	EA18,
	EA19,
	EA20,
	EA21,
	EA22,
	EA23,
	EA24,
	EA25,
	EA26,
	EA27,
	EA28,
	EA29,
	EA30,
	EA31,
	EA32,
	EA33,
	EA34,
	EA35,
	EA36,
	EA37,
	EA38,
	EA39,
	EA40,
	EA41,
	EA42,
	EA43,
	EA44,
	EA45,
	EA46,
	EA47,
	EA48,
	EA49,
	EA50,
	EA51,
	EA52,
	EA53,
	EA54,
	EA55,
	EA56,
	EEndNuclei,
	EEndParticle = EEndNuclei,
	EBeginReservedParticles = EEndParticle,

	//reserved for future (to preserve output spectrum file format)
	EReserved1 = EBeginReservedParticles,
	EReserved2,
	EReserved3,
	EReserved4,
	EReserved5,
	EReserved6,
	EEndReservedParticles,
	
	EEndAllParticles = EEndReservedParticles,
	
	ENucleon=ENeutron,
	EEndEM=ENeutrinoE//end electromagnetic
};

#ifdef _DEBUG_TPARTICLE
class TParticle{
public:
	TParticle(int aParticle){
		m_particle = (_TParticle)aParticle;
		int check = aParticle - (int)m_particle;
		ASSERT(check==0);
	}
	operator int() const {return (int)m_particle;}
private:
	_TParticle m_particle;
};
#else // #ifdef _DEBUG_TPARTICLE
typedef _TParticle TParticle;
#endif

#endif //#if !defined(TPARTICLE_H_INCLUDED)
