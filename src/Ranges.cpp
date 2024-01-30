// Ranges.cpp: implementation of the CRanges class.
//
//////////////////////////////////////////////////////////////////////

#include "Ranges.h"
#include <math.h>
#include "ParticleData.h"
#include "Units.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//double Emax=1.0e17;    /*Max. input energy of ray (MeV) */
//double Emin=1.0e2;    /*Min. final energy of ray (MeV) */
//double emin=0;    /*Current min. energy of ray (Eunits) */
//double s=20.0; /* average number of points per decade */ 
//int nn=56;   /* number of points between Emin & Emax*/
//int n10 = 0; /* if not equal to 0 defines number of points per decade, and thus defines nn*/


//double ss1,ss2,ss_1,ss_2,ss3,sF;  /* constants which depend on s */

const int CRanges::s_minNumberOfBinsPerDecade=10;

CRanges::CRanges(int aAccuracy, double aApproxEmin, double aApproxEmax, 
		//double aApproxKmin, double aApproxKmax,
		double aZmax, int aAccuracyZ, bool aCELredshift, double aMinInteractionEnergyMeV,
		double aMinNucleiInteractionEnergyMeV):
m_Emin(aApproxEmin),
m_Emax(aApproxEmax),
m_bc(aAccuracy*s_minNumberOfBinsPerDecade),
m_Zmax(aZmax),
m_accuracy(aAccuracy),
m_accuracyZ(aAccuracyZ)
{
	if(m_Emin>=m_Emax)
		ThrowError("Invalid range (Emin>=Emax)");
	m_nE = adjustRanges(m_Emin, m_Emax);
	//m_nK = adjustRanges(m_Kmin, m_Kmax);

	m_energies.Create(m_Emin, m_nE, m_bc);
	m_midEnergies.Create(m_Emin*m_bc.ss2, m_nE, m_bc);
	m_nucleonGammas.Create(m_Emin/ParticleData::meanNucleonMassMeV*units.Eunit, m_nE, m_bc);
	m_midNucleonGammas.Create(m_Emin*m_bc.ss2/ParticleData::meanNucleonMassMeV*units.Eunit, m_nE, m_bc);
	//m_backgroundEnergies.Create(m_Kmin, m_nK, m_bc);
	//m_midBackgroundEnergies.Create(m_Kmin*m_bc.ss2, m_nK, m_bc);
	m_MinInteractionBin = m_energies.getBin(aMinInteractionEnergyMeV/units.Eunit);
    m_MinInteractionBinN = m_energies.getBin(aMinNucleiInteractionEnergyMeV/units.Eunit);
	if(m_MinInteractionBin < 0)
		m_MinInteractionBin = 0;
	if(m_MinInteractionBinN < m_MinInteractionBin)
        m_MinInteractionBinN = m_MinInteractionBin;

	//initialize z - binning
	{
		double coef = pow(m_bc.ss1, 1./((double)m_accuracyZ));//m_accuracyZ-th part of energy bin
		if(!aCELredshift)
		{
			int nBins = m_accuracyZ * ((int)(log(1.+m_Zmax)/log(m_bc.ss1)+1.));
			m_Z.create(nBins+1);
			double z = 0.;
			for(int i=0; i<=nBins; i++)
			{
				m_Z[i] = z;
				z=coef*(1.+z)-1.;
			}
			m_Zmax = m_Z[nBins];
		}
		else
		{
			int nBins = 0;
			double z = 0.;
			for(;z<m_Zmax;nBins++)
				z=coef*(1.+z)-1.;
			m_Z.create(nBins+1);
                        std::cerr << "Z-binning:\n";
			for(z=0.0,nBins=0;z<m_Zmax;nBins++)
			{
				m_Z[nBins] = z;
				z=coef*(1.+z)-1.;
                                std::cerr << nBins << "\t" << m_Z[nBins] << "\n";
			}
			m_Z[nBins] = m_Zmax;
                        std::cerr << nBins << "\t" << m_Z[nBins] << std::endl;
		}
                
	}
}

/// creating binning with exactly iAccuracy * s_minNumberOfBinsPerDecade points per decade 
/// in the way to ensure that for any integer n E=10^n eV fits exactly into the middle of a bin
int CRanges::adjustRanges(double& aMin, double& aMax) const{

	aMin *= units.Eunit;// converting to MeV
	aMax *= units.Eunit;// converting to MeV

	//double noOfBinsPerDecade = m_bc.s;

	int n = (int)(m_bc.s*log10(aMin)+0.5);
	aMin = pow(10.,((double)n - 0.5)/m_bc.s);

	int m = (int)(m_bc.s*log10(aMax)+1.5);
	aMax = pow(10.,((double)m - 0.5)/m_bc.s);

	aMin /= units.Eunit;// converting back to internal units
	aMax /= units.Eunit;// converting back to internal units
	
	return m-n;
}

CRanges::~CRanges()
{

}

CRanges*  CRanges::s_defaultRanges = NULL;

CRanges* CRanges::SetDefault(CRanges* aRanges){
	CRanges*  result = s_defaultRanges;
	s_defaultRanges = aRanges;
	return result;
}

TBinningConstants::TBinningConstants(int aNoOfBinsPerDecade){
	NoOfBinsPerDecade = aNoOfBinsPerDecade;
	s = aNoOfBinsPerDecade;
	ss1=pow(10.0,1.0/s);
	ss2=sqrt(ss1);
	ss_1=1.0/ss1;
	ss_2=1.0/ss2;
	ss3=ss1*ss1*ss1;
	sF=ss2-ss_2;
	sLn = log(ss1);
}

void CBinning::Create(double aInitialVal, int aNoOfBins, const TBinningConstants& aBc)
{
	create(aNoOfBins);
	double x = aInitialVal;
	for(int i=0;i<aNoOfBins;i++,x*=aBc.ss1)
	{
		(*this)[i] = x;
	}
	m_minVal = (*this)[0]/aBc.ss2;
	m_maxVal = (*this)[m_length-1]*aBc.ss2;
}

