// Vector.cpp: implementation of the CVector class.
//
//////////////////////////////////////////////////////////////////////

#include <istream>
#include <iterator>
#include "Vector.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

CVector& CVector::operator=(const CBaseVector& theVector)
{
//	ASSERT(length()==theVector.length());
//	ASSERT(getIndexShift()==theVector.getIndexShift());
//	for(int i=m_iMin;i<=m_iMax;i++)
//		(*this)[i] = theVector[i];
	CBaseVector::copy(theVector);
	return (*this);
}

CVector& CVector::operator+=(const CBaseVector& theVector)
{
	ASSERT(length()==theVector.length());
	ASSERT(getIndexShift()==theVector.getIndexShift());
	for(int i=m_iMin;i<=m_iMax;i++)
		(*this)[i] += theVector[i];
	return (*this);
}

CVector& CVector::operator-=(const CBaseVector& theVector)
{
	ASSERT(length()==theVector.length());
	ASSERT(getIndexShift()==theVector.getIndexShift());
	for(int i=m_iMin;i<=m_iMax;i++)
		(*this)[i] -= theVector[i];
	return (*this);
}

CVector& CVector::operator*=(double theMult)
{
	for(int i=m_iMin;i<=m_iMax;i++)
		(*this)[i] *= theMult;

	return (*this);
}

void CVector::reset()
{
	for(int i=m_iMin;i<=m_iMax;i++)
		(*this)[i] = 0.;
}

CVector& CVector::operator/=(double theMult)
{
	ASSERT(theMult!=0);
	double mult = 1./theMult;
	(*this) *= mult;
	return (*this);
}

int CVector::findX(double xValue, double* xArray, int xArraySize)
//find the nearest to xValue x-array element (increasing x-array supposed)
//if value is out of range, returns -1
//using dihotomy method
{
	if((xValue>xArray[xArraySize-1])||xValue<xArray[0])
		return -1;
	int left = 0;
	int right = xArraySize-1;
	for(int i=(left+right)/2;right-left>1;i=(left+right)/2)
	{
		if(xValue>xArray[i])
			left=i;
		else
			right=i;
	}
	ASSERT((right - left) == 1);
	return ((xArray[right]-xValue)<(xValue-xArray[left]))?right:left;
}

//void CVector::copy(const double* v)
//{
//	ASSERT(m_iMax>m_iMin);//array must not be empty
//	//memcpy(m_array+m_iMin,v,sizeof(double))
//todo
//	int i,j;
//	for(i=m_iMin,j=0;i<=m_iMax;i++,j++)
//		(*this)[i] = v[j];
//}

int CVector::findLeftX(double xValue,double* xArray,int xArraySize)
{
	if((xValue>xArray[xArraySize-1])||xValue<xArray[0])
		return -1;
	int left = 0;
	int right = xArraySize-1;
	for(int i=(left+right)/2;right-left>1;i=(left+right)/2)
	{
		if(xValue>xArray[i])
			left=i;
		else
			right=i;
	}
	ASSERT((right - left) == 1);
	return left;
}

#ifdef USE_GSL
gsl_vector* CVector::GslVector() const
{
	gsl_vector* x = gsl_vector_alloc (m_length);
	for(int i=0; i<m_length; i++)
		gsl_vector_set (x, i, m_array[i]);
	return x;
}

void  CVector::create_gsl(const gsl_vector *v)
{
	create(v->size);
	for(int i=v->size - 1; i>=0; i--)
		m_array[i] = gsl_vector_get(v, i);
}

void  CVector::copy_gsl(const gsl_vector *v)
{
	ASSERT(m_length == (int)v->size);
	for(int i=m_length - 1; i>=0; i--)
		m_array[i+m_iMin] = gsl_vector_get(v, i);
}
#endif

//find the nearest to xValue x-array element (increasing x-array supposed)
//if value is out of range, returns <minimal array index>-1
int CVector::findX(double xValue) const
{
	return findX(xValue, m_array + m_iMin, m_length) + m_iMin;
}

//find the nearest to xValue x-array element which is less than the argument (increasing x-array supposed)
int CVector::findLeftX(double xValue) const
{
	return findLeftX(xValue, m_array + m_iMin, m_length) + m_iMin;
}

CMatrix::CMatrix()
{

}

CMatrix::CMatrix(std::istream& aInput)
{
	std::string line="";
	int nCols = 0;
	CAutoDeletePtrArray<CVarVector>  data;
	for(std::getline(aInput, line);!aInput.eof();std::getline(aInput, line))
	{
		trim(line);
		if(line.length()==0 || line[0]=='#')
			continue;
		std::istringstream input(line);
		std::istream_iterator<double> eos;
		std::istream_iterator<double> iit (input);
		CVarVector* row = new CVarVector();
		data.add(row);

		int recNo=0;

		for(; iit!=eos; recNo++, iit++)
		{
			row->add(*iit);
		}
		if(nCols && nCols!=recNo)
			ThrowError("Matrix reading error - variable number of columns");
		nCols = recNo;
	}

	int iSize = data.length();
	if(iSize==0)
		return;
	m_array.create(iSize);
	for(int i=0;i<iSize;i++)
		m_array[i] = new CVector(data(i));
	m_type = (iSize == data[0]->length())?squareE:rectangleE;
}

void CMatrix::create(int iSize, int jSize)
{
	m_array.create(iSize);
	for(int i=0;i<iSize;i++)
		m_array[i] = new CVector(jSize);
	m_type = (iSize==jSize)?squareE:rectangleE;
}

void CMatrix::multiply(double aValue)
{
	int iMin = m_array.getIndexShift();
	int iMax = m_array.length()+iMin;
	for(int i=iMin; i<iMax; i++)
		(*(m_array[i])) *= aValue;
}

void CMatrix::solveLinearSystem (CVector& b) const
{
	int size = m_array.length();
	ASSERT(size == b.length());
	gsl_matrix* A = gsl_matrix_alloc (size, size);
	gsl_vector* x = gsl_vector_alloc(size);
	for(int i=0; i<size; i++)
	{
		gsl_vector_set(x, i, b[i]);
		for(int j=0; j<size; j++)
			gsl_matrix_set(A, i, j, (*this)[i][j]);
	}
	int error = gsl_linalg_HH_svx (A, x);
	if(error)
	{
		ASSERT(0);
		gsl_vector_free(x);
		gsl_matrix_free(A);
		ThrowError("Failed to solve linear system: error code " + ToString(error));
	}
	for(int i=0; i<size; i++)
		b[i] = gsl_vector_get(x, i);
	gsl_vector_free(x);
	gsl_matrix_free(A);
}

void CMatrix::reset()
{
	int iMin = m_array.getIndexShift();
	int iMax = iMin + m_array.length();
	for(int i=iMin; i<iMax; i++)
		m_array[i]->reset();
}

void CMatrix::create(int size, CMatrix::EMatrixType type)
{
	m_array.create(size);
	switch(type)
	{
	case squareE:
	case rectangleE:
		{
			for(int i=0;i<size;i++)
				m_array[i] = new CVector(size);
				type = squareE;
		}
		break;
	case triangleLeftE:
		{
			for(int i=0;i<size;i++)
				m_array[i] = new CVector(i+1);
		}
		break;
	case triangleRightE:
		{
			for(int i=0;i<size;i++)
			{
				m_array[i] = new CVector(size-i);
				m_array[i]->setIndexShift(-i);
			}
		}
		break;
	}
	m_type = type;
}

CMatrix::CMatrix(int iSize, int jSize)
{
	create(iSize,jSize);
}

CMatrix::CMatrix(int size, CMatrix::EMatrixType type)
{
	create(size,type);
}

CMatrix::~CMatrix()
{
	int i = - m_array.getIndexShift();
	int iMax = i + m_array.length();
	for(; i<iMax ;i++)
		delete m_array[i];
}
