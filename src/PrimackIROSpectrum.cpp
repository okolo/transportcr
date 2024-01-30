#include "PrimackIROSpectrum.h"
#include <math.h>

#define IRO_DATA_FILE "iro.dat"

CPrimackIROSpectrum::CPrimackIROSpectrum():
CTableWithHeaderBackgr(DATA_DIR IRO_DATA_FILE),
m_zAc(0),
m_eAc(0),
m_eBuffer(MAX_IRO_Z_ACCURASSY),
m_nBuffer(MAX_IRO_Z_ACCURASSY)
{
	
}

CPrimackIROSpectrum::~CPrimackIROSpectrum()
{

}

CPrimackIROSpectrum::ECompleteStatus CPrimackIROSpectrum::readHeaderLine(const char* theString)
{


	if(m_curHeaderLine==1)
	{
		if(sscanf(theString," zAc = %d eAc = %d",&m_zAc, &m_eAc)!=2)
			return failedE;
		if((m_zAc<2)||(m_eAc<2))
			return failedE;
		ASSERT(m_zAc<=MAX_IRO_Z_ACCURASSY);
		m_zArray.create(m_zAc);
		m_eArray.create(m_zAc,m_eAc);
		m_nArray.create(m_zAc,m_eAc);
		m_eAlphaArray.create(m_zAc,m_eAc);
		return notCompletedE;
	}
	else if(m_curHeaderLine==2)
	{
		if(parseZvaluesString(theString))
			return completedE;
		else
			return failedE;
	}
	ASSERT(false);
	return failedE;
}

bool CPrimackIROSpectrum::readDataLine(const char* theString)
{
	ASSERT(m_zAc<=MAX_IRO_Z_ACCURASSY);
	char formatStr[512];
	formatStr[0]='\0';

#define ENTRY(i) (m_eBuffer + i),(m_nBuffer + i)
	int i;
	for(i=0;i<m_zAc;i++)
		strcat(formatStr," %lg %lg");

	// here it's important that m_zAc < MAX_IRO_Z_ACCURASSY !!!
	if(sscanf(theString,formatStr,ENTRY(0),ENTRY(1),ENTRY(2),
		ENTRY(3),ENTRY(4),ENTRY(5))!=2*m_zAc)
		return false;
	for(i=0;i<m_zAc;i++)
	{
		m_eArray[i][m_curDataLine-1] = m_eBuffer[i];
		m_nArray[i][m_curDataLine-1] = m_nBuffer[i];
	}
#undef ENTRY
	return true;
}

bool CPrimackIROSpectrum::testData()
{
	return true;
}

void CPrimackIROSpectrum::processData()
{
	for(int i=0;i<m_zAc;i++)
	{
		if(m_eArray[i][0]>m_eArray[i][m_eAc-1])
		{
			m_eArray[i].invert();
			m_nArray[i].invert();
		}
		for(int j=0;j<m_eAc;j++)
		{
			if(j<m_eAc-1)
				m_eAlphaArray[i][j]=(m_nArray[i][j+1]-m_nArray[i][j])/(m_eArray[i][j+1]-m_eArray[i][j])-1.;
			double x = pow(10.,m_eArray[i][j]);
			double y = pow(10.,m_nArray[i][j])/x;
			m_eArray[i][j] = x;
			m_nArray[i][j] = y;
		}
	}
}


bool CPrimackIROSpectrum::parseZvaluesString(const char *theString)
{
	ASSERT(m_zAc<=MAX_IRO_Z_ACCURASSY);
#define BUF(i) (m_eBuffer + i)
	char formatStr[256];
	formatStr[0]='\0';
	int i;
	for(i=0;i<m_zAc;i++)
		strcat(formatStr," z = %lg");
	if(sscanf(theString,formatStr,BUF(0),BUF(1),BUF(2),
		BUF(3),BUF(4),BUF(5))!=m_zAc)
		return false;
	for(i=0;i<m_zAc;i++)
		m_zArray[i] = m_eBuffer[i];
	ASSERT(m_zArray[0]==0.);
	return true;
#undef BUF
}

double CPrimackIROSpectrum::MinE(double aZmax) const
{
	return m_eArray[0][0];
}

double CPrimackIROSpectrum::MaxE(double aZmax) const
{
	return m_eArray[0][m_eAc-1];
}

double CPrimackIROSpectrum::MaxZ() const
{
	return m_zArray[m_zAc-1];
}

double CPrimackIROSpectrum::F(double E, double z)
{
	UINT iZ=0;
	int iE=-1;
	if(z<1e-4)//z=0
	{
		return F(E,iZ);
	}
	int intZ=m_zArray.findLeftX(z);
	if(intZ<0)
		return 0.;
	iZ = (UINT)intZ;
	double left = F(E,iZ,&iE);
	if(iE<0)
		return 0.;

	double xlog = log((m_zArray[iZ+1]+1.)/(m_zArray[iZ]+1.));
	double right = F(E,iZ+1);
	if(right<=0.)
		return 0.;
	double alpha = log(right/left)/xlog;//alpha is power law parameter in (1+z) dependance of F
	return left*pow((z+1.)/(m_zArray[iZ]+1.),alpha);
		// + (right-left)/(m_zArray[iZ+1]-m_zArray[iZ])*(z-m_zArray[iZ]);
}

double CPrimackIROSpectrum::F(double E, UINT iZ, int* iE)
{
	int i=m_eArray[iZ].findLeftX(E);
	if(iE!=NULL)
		*iE = i;
	if(i<0)
		return 0.;
	return m_nArray[iZ][i]*pow(E/m_eArray[iZ][i],m_eAlphaArray[iZ][i]);
}

double CPrimackIROSpectrum::F(double E, UINT iZ)
{
	if ((int)iZ>=m_zAc) ThrowError("CPrimackIROSpectrum::F() : invalid argument");
	return F(E, iZ, NULL);
}

