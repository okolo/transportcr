// Stecker2005IROSpectrum.cpp: implementation of the CStecker2005IROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "Stecker2005IROSpectrum.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#define IRO_DATA_FILE "iro_stecker2005"
CStecker2005IROSpectrum::CStecker2005IROSpectrum():
CTableWithHeaderBackgr(DATA_DIR IRO_DATA_FILE),
m_accuracyE(0),
m_accuracyZ(0),
m_curE(0)
{

}

CStecker2005IROSpectrum::~CStecker2005IROSpectrum()
{

}

double CStecker2005IROSpectrum::F(double E, double z)
{
	int iZ=0;
	//	int iE=-1;
	if(z<1e-4)//z=0
	{
		return m_fE[0]->f(E);// logarifmic approximation on E
	}
	iZ=m_zArray.findLeftX(z);
	if(iZ<0)
		return 0.;

	double n1 = m_fE[iZ]->f(E);// logarifmic approximation on E
	double n2 = m_fE[iZ+1]->f(E);
	double z1 = m_zArray[iZ];
	double z2 = m_zArray[iZ+1];

	double result = n1 +  (n2-n1)/(z2-z1)*(z-z1);  // simple linear approximation on z

	return result;
}

double CStecker2005IROSpectrum::MinE(double aZmax) const
{
	return m_eArray[0];
}

double CStecker2005IROSpectrum::MaxE(double aZmax) const
{
	return m_eArray[m_accuracyE-1];
}

double CStecker2005IROSpectrum::MaxZ() const
{
	return m_zArray[m_accuracyZ-1];
}

CTableWithHeaderReader::ECompleteStatus CStecker2005IROSpectrum::readHeaderLine(const char* theString)
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

bool CStecker2005IROSpectrum::readDataLine(const char* theString)
{
	if (m_curE>=m_accuracyE)
	{
		return false;
	}
	CVector v;
	if(!readVector(v, theString))
		return false;
	ASSERT(v.length()==m_accuracyZ+1);
	double E = pow(10.,v[0]-9.)/241803.; // first column contains log10(E/Hz), converting to eV
	m_eArray[m_curE] = E;
	for(int i=0; i< m_accuracyZ; i++)
		//m_nArray[i][m_curE] = v[i+1]*1e-6*E*pow(1.+m_zArray[i],-3); // multiply by energy and convert m^-3 to cm^-3, multiplied by (1+z)^{-3} to make it in comoving volume
		m_nArray[i][m_curE] = v[i+1]*1e-6*pow(1.+m_zArray[i],-3); // convert m^-3 to cm^-3, multiplied by (1+z)^{-3} to make it in comoving volume
	m_curE++;
	return true;
}

bool CStecker2005IROSpectrum::testData()
{
	return true;
}

void CStecker2005IROSpectrum::processData()
{
	for(int i=0; i<m_accuracyZ; i++)
		m_fE[i] = new CDefaultTableFunc(m_eArray, m_nArray[i]);
}

Stecker2012IROSpectrum::Stecker2012IROSpectrum(bool aUpperLevel):
		f2012specUV(aUpperLevel?"stecker2012/upperdataIV.dat":"stecker2012/lowerdataIV.dat",
				false, false),
				fMaxZ(0)
{}

bool Stecker2012IROSpectrum::init()
{
	bool result = f2005spec.init();
	result = f2012specUV.init() && result;
	fMaxZ = f2012specUV.MaxZ() > f2005spec.MaxZ() ? f2005spec.MaxZ() : f2012specUV.MaxZ();
	return result;
}

/* IR/O spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
	  (must be multiplied by (1+z)^3 before substituting farther)*/
double Stecker2012IROSpectrum::F(double E, double z)
{
	return E>1.58489319246111 ? f2012specUV.F(E,z) : f2005spec.F(E,z);
}

double Stecker2012IROSpectrum::MaxZ() const
{
	return fMaxZ;
}

double Stecker2012IROSpectrum::MaxE(double aZmax) const
{
	return f2012specUV.MaxE(aZmax);
}

double Stecker2012IROSpectrum::MinE(double aZmax) const
{
	return f2005spec.MinE(aZmax);
}
