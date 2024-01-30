// KneiskeIROSpectrum.cpp: implementation of the CKneiskeIROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "KneiskeIROSpectrum.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#define IRO_DATA_FILE "EBL_kneiske.dat"
CKneiskeIROSpectrum::CKneiskeIROSpectrum():
CTableWithHeaderBackgr(DATA_DIR IRO_DATA_FILE),
m_accuracyE(0),
m_accuracyZ(0),
m_curE(0)
{

}

CKneiskeIROSpectrum::~CKneiskeIROSpectrum()
{

}

double CKneiskeIROSpectrum::F(double E, double z)
{
	int iZ=0;
	//	int iE=-1;
	if(z<1e-4)//z=0
	{
		return m_fE[0]->f(E);// logarithmic approximation on E
	}
	iZ=m_zArray.findLeftX(z);
	if(iZ<0)
		return 0.;

	double n1 = m_fE[iZ]->f(E);// logarithmic approximation on E
	double n2 = m_fE[iZ+1]->f(E);
	double z1 = m_zArray[iZ];
	double z2 = m_zArray[iZ+1];

	double result = n1 +  (n2-n1)/(z2-z1)*(z-z1);  // simple linear approximation on z

	return result;
}

double CKneiskeIROSpectrum::MinE(double aZmax) const
{
	return m_eArray[0];
}

double CKneiskeIROSpectrum::MaxE(double aZmax) const
{
	return m_eArray[m_accuracyE-1];
}

double CKneiskeIROSpectrum::MaxZ() const
{
	return m_zArray[m_accuracyZ-1];
}

CTableWithHeaderReader::ECompleteStatus CKneiskeIROSpectrum::readHeaderLine(const char* theString)
{
	if (m_accuracyE<=0)
	{
		return readDataLength(m_accuracyE,theString)?notCompletedE:failedE;
	}
	
	if (!readVector(m_zArray, theString))
	{
		return failedE;
	}

	m_eArray.create(m_accuracyE);
	m_accuracyZ = m_zArray.length();
	m_nArray.create(m_accuracyZ, m_accuracyE);
	m_fE.create(m_accuracyZ);

	return completedE;
}

bool CKneiskeIROSpectrum::readDataLine(const char* theString)
{
	if (m_curE>=m_accuracyE)
	{
		return false;
	}
	CVector v;
	if(!readVector(v, theString))
		return false;
	ASSERT(v.length()==m_accuracyZ+1);
	double E = (m_eArray[m_curE] = v[0]); // first column contains value of energy
	for(int i=0; i< m_accuracyZ; i++)
		m_nArray[i][m_curE] = v[i+1]*E*1e-6; // multiply by energy and convert m^-3 to cm^-3

	m_curE++;
	return true;
}

bool CKneiskeIROSpectrum::testData()
{
	return true;
}

void CKneiskeIROSpectrum::processData()
{
	for(int i=0; i<m_accuracyZ; i++)
		m_fE[i] = new CDefaultTableFunc(m_eArray, m_nArray[i]);
}
