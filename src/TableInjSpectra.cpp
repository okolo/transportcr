// UserInjSpectra.cpp: implementation of the CTableInjSpectra class.
//
//////////////////////////////////////////////////////////////////////

#include "TableInjSpectra.h"
#include "main.h"
//#include "Prodspec.h"
#include "ParticleList.h"
#include "Concentrations.h"
#include "Parameters.h"
#include "Addfunc.h"
#include "Units.h"
#include "PropagEngine.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTableInjSpectra::CTableInjSpectra()
{
	m_DataFileName = "all";
	Reader()->readStringPar("CustomInjectionFileName", m_DataFileName);
	if(CPropagEngine::GetSourceType() == PointSource) {
		m_deltaT = SourceTable::DefaultMomentSourceDt;//make sure initial spectrum has same norm as input file is absolute norm is used and SourceAbsNorm=1
	}
	else {
		double CustomInjectionDeltaT_Mpc = 1.;//Mpc
		READ_DOUBLE_SETTING(CustomInjectionDeltaT_Mpc);
		m_deltaT = CustomInjectionDeltaT_Mpc * units.Mpc;
	}
}

void CTableInjSpectra::init()
{
	m_spectra = new Concentrations();
	std::string path = "custom_injection";
	path = path + DIR_DELIMITER_STR + m_DataFileName;
	CopyFileOrDir(path ,plt_local_dir + m_DataFileName);
	m_spectra->ReadAll(path);
}

double CTableInjSpectra::Q(TParticle aParticle, int aBinE, double aEinMeV, double aZ)
{
	double coef=aEinMeV*BC().sF*m_deltaT/(units.sec*units.cm3);
	ASSERT(coef>0 && coef < 1e300);
	//we assume that concentration in m_spectra corresponds to source emitting for time m_deltaT
	if(aBinE>=0)
		return m_spectra->Bin(aParticle, aBinE)/coef;
	else
		return m_spectra->GetValue(aParticle, aEinMeV * units.MeV)/coef;
}

