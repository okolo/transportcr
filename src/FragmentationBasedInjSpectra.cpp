#include "FragmentationBasedInjSpectra.h"
#include "TableReader.h"
#include "Parameters.h"
#include "ParticleList.h"
#include "TableFunc.h"

CFragmentationBasedInjSpectra::CFragmentationBasedInjSpectra():
m_graphPower(3),
m_SourceDir("custom_fragmentation")
{
	Reader()->readDoublePar("FragmentationGraphPower", m_graphPower);
}

CFragmentationBasedInjSpectra::~CFragmentationBasedInjSpectra(void)
{
}

double CFragmentationBasedInjSpectra::Q(TParticle aParticle, int aBinE, double aE/*MeV*/, double aZ)
{
	double maxE = MaxE(EProton);
	CTableFunction* fragmentationGraph = m_fragSpectra[aParticle];

	if(fragmentationGraph == NULL)
		return 0.;
	
	double x = aE/maxE;

	double result = fragmentationGraph->fApprox(x);/// also approximating values beyond table assuming linearity in log scale 

	if(result<=0.) return 0.;
	
	result *= pow(x, -m_graphPower);

//	weight*BC().sF*E*SpecUnit*Q##particle##0(Eunit*E,i) - note, farther the code like this is executed
	return result/maxE;
}

void CFragmentationBasedInjSpectra::init()
{
	m_fragSpectra.create(EEndParticle);
	FOR_ALL_REAL_PARTICLES_INVOLVED(particle)
	{
		const char* fileName = ParticleData::getParticleFileName(particle);
		FILE* file = Fopen(fileName, m_SourceDir, "rt", false);
		m_fragSpectra[particle] = file?(new PropagDepreciated::CLogScaleFunc(new CTableReader(file,2))):NULL;
	}
}

