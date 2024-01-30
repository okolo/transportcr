// LEPData.cpp: implementation of the CLEPData class.
//
//////////////////////////////////////////////////////////////////////

#include "LEPData.h"
#include "Addfunc.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#define LEP_DATA_FILE "LEPdata.txt"

CLEPData::CLEPData():
CTableWithHeaderReader(DATA_DIR LEP_DATA_FILE),
x(NULL),
y(NULL),
m_accuracy(0),
m_noOfDecays(0),
m_currentPoint(0)
{
}

CLEPData::~CLEPData()
{
	delete[] x;
	if(y!=NULL)
		for(int i=0;i<(int)endE;i++)
			delete[] (y[i]);
	delete[] y;
}

CLEPData::ECompleteStatus CLEPData::readHeaderLine(const char* theString)
{
	if(sscanf(theString,
		" accuracy = %u , total = %u",
		&m_accuracy,&m_noOfDecays)==2)
	{
		x = new double[m_accuracy + 1];
		x[m_accuracy]=1.;
		y = new PDouble[(int)endE];
		for(int i=0;i<(int)endE;i++)
		{
			y[i] = new double[m_accuracy+1];
			y[i][m_accuracy] = 0;
		}

		m_currentPoint = 0;
		m_controlSum = 0.;
		return completedE;
	}
	else
		return failedE;
}

bool CLEPData::readDataLine(const char* theString)
{
	if(m_currentPoint>=m_accuracy)
		return false;
	double right;
	int curIndex = m_accuracy - 1 - m_currentPoint;
	int i = sscanf(theString," %lg-%lg %lg %lg %lg %lg %lg"
		,&right
		,x + curIndex
		,y[electronsE] + curIndex
		,y[photonsE] + curIndex
		,y[protonsE] + curIndex
		,y[neutronsE] + curIndex
		,y[neutrinosE] + curIndex);
//energy conservation test
	double E = pow(10.,-0.5*(x[curIndex]+right));// 2*E/m_z
	for(int j=0;j<endE;j++)
		m_controlSum+=(y[j][curIndex])*E;

	m_currentPoint++;
	return ((i-2)==endE);
}

bool CLEPData::testData()
{
	m_controlSum/=(2*m_noOfDecays);//must be <~1
//	ASSERT(m_controlSum<1.);//==Gamma_hadron/Gamma_tot
//	printf("\nLEP sum = %lg\n",m_controlSum);
	return (m_currentPoint==m_accuracy);
}

void CLEPData::processData()
{
  //double ln10 = log(10.);
	double xPrev = 1.;
	for(int i=m_accuracy-1;i>=0;i--)
	{
		x[i] = pow(10.,-x[i]);
		double xMiddle= sqrt(x[i]*xPrev);
		xPrev = x[i];
		for(int particle = 0 ; particle<endE ; particle++ )
		{
			y[particle][i] /= xMiddle*m_noOfDecays;
			if(i<((int)(m_accuracy))-1)
				y[particle][i] += y[particle][i+1];
		}
	}
	DEBUG_ONLY(Save());
}

double CLEPData::N(EParticles particle,double _x)
{
	ASSERT(particle!=endE);

	double result1 = yValue(_x, x, y[particle], m_accuracy + 1, y[particle][0], 0.);
/*	double result2 = (_x<1.)?(m_controlSum*(1./_x-1.)/5.):0.;//test
	if(result2<result1)
		return result1;*/
	return result1;
}

void CLEPData::Save()
{
	FILE *file = Fopen("CLEPData.N",plt_local_c,"wt",true);
	fprintf(file,"#2E/Mz      #electrons    #photons    #protons   #neutrons  #neutrinos");
	for(int i=0;i<(int)m_accuracy;i++)
	{
		fprintf(file,"\n%lg",x[i]);
		for(int particle=0;particle<endE;particle++)
			fprintf(file,"  %lg",y[particle][i]);
	}
	fclose(file);
}
