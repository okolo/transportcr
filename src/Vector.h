// Vector.h: interface for the CVector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
#define AFX_VECTOR_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Addfunc.h"
#include <vector>
#ifdef USE_GSL
#include <gsl/gsl_vector_double.h>
#endif

template<class type, class argType>
class Array  : public TReferencedObj
{
public:
	Array(int size=0):
	m_length(0),
	m_array(NULL),
	m_iMin(0),
	m_iMax(-1),//m_iMin>m_iMax indicates that array is empty
	m_OwnsMemory(true)
	{
		ASSERT(size>=0);
		if(size>0)
			create(size);
	};

	Array(const std::vector<type>& std_vec):Array(std_vec.size()){
        memcpy(m_array, std_vec.data(),m_length*sizeof(type));
	}

	Array(const Array<type,argType>& array):
	m_length(array.m_length),
	m_array(NULL),
	m_iMin(array.m_iMin),
	m_iMax(array.m_iMax),
	m_OwnsMemory(true)
	{
		if(array.length()<=0)
			return;
		m_array = (type*)Calloc(m_length,sizeof(type),"Array");
		memcpy(m_array,array.m_array + m_iMin,m_length*sizeof(type));
		m_array-=m_iMin;
	}

	Array(type* aArray, int aLength):
	m_length(aLength),
	m_array(aArray),
	m_iMin(0),
	m_iMax(aLength-1),
	m_OwnsMemory(false)
	{

	}

	void appendToSTLVector(std::vector<type>& aSTLvector)
	{
		for(int i=m_iMin; i<=m_iMax;i++)
			aSTLvector.push_back(m_array[i]);
	}

	void create(const Array<type,argType>& array)
	{
		if(array.length()<=0)
			return;
		ASSERT(m_length==0);//create() method can be called only once
		m_length = array.length();
		m_iMin = array.m_iMin;
		m_iMax = array.m_iMax;
		m_array = (type*)Calloc(m_length,sizeof(type),"Array");
		m_OwnsMemory = true;
		memcpy(m_array,array.m_array + m_iMin,m_length*sizeof(type));
		m_array-=m_iMin;	
	};
		
	virtual void create(int size)
	{
		ASSERT(size>0);
		ASSERT(m_length==0);//create() method can be called only once
		m_length = size;
		m_iMax = size-1;
		m_array = (type*)Calloc(size,sizeof(type),"Array");
		m_OwnsMemory = true;
	};

	virtual ~Array()
	{
		if(m_OwnsMemory)
			free(m_array+m_iMin);
	};

	virtual void copy(const Array<type, argType>& aFrom){//binary copy of elements
		ASSERT(aFrom.length()==length());
		ASSERT(aFrom.getIndexShift()==getIndexShift());
		memcpy(m_array+m_iMin, aFrom.m_array+m_iMin,sizeof(type)*m_length);
	}

	virtual void copy(const type* aFrom){//binary copy of elements
		ASSERT(m_length > 0);
		memcpy(m_array+m_iMin, aFrom, sizeof(type)*m_length);
	}

	Array<type, argType>& operator = (const Array<type, argType>& aFrom){
		ASSERT(aFrom.length()==length());
		ASSERT(aFrom.getIndexShift()==getIndexShift());
		for(int i=m_iMin; i<=m_iMax;i++)
			m_array[i] = aFrom[i];
		return (*this);
	}
	
	void invert()
	{
		type* pArray = m_array+m_iMin;
		int maxIndex = m_length/2;
		int i1,i2;
		for(i1=0,i2=m_length-1;i1<maxIndex;i1++,i2--)
		{
			double mem = pArray[i1];
			pArray[i1] = pArray[i2];
			pArray[i2] = mem;
		}
	};
	
	virtual void moveIndex(int theShift)
	{
		ASSERT(m_length>0 || theShift==0);
		m_array+=theShift;
		m_iMin-=theShift;
		m_iMax-=theShift;
	};
		
	virtual void setIndexShift(int theShift)
	{
		ASSERT(m_length>0 || theShift==0);
		
		// moving to shift 0
		m_array += m_iMin;
		
		// adjusting shift
		m_iMax = m_length - 1 - theShift;
		m_iMin = -theShift;
		
		// moving to new shift
		m_array -= m_iMin;		
	};

	Array<type,argType>& operator+=(const int i)
	{
		moveIndex(i);
		return *this;
	};

	Array<type,argType>& operator++()
	{
		moveIndex(1);
		return *this;
	}

	Array<type,argType>& operator-=(const int i)
	{
		moveIndex(-i);
		return *this;
	};

	Array<type,argType>& operator--()
	{
		moveIndex(-1);
		return *this;
	}

	inline int getIndexShift() const {return -m_iMin;};

	inline argType operator[](int i)
	{
		ASSERT((i>=m_iMin)&&(i<=m_iMax));
		return m_array[i];
	};
	
	inline const argType operator[](int i) const
	{
		ASSERT((i>=m_iMin)&&(i<=m_iMax));
		return m_array[i];
	};
	
	inline type* operator+(int i)
	{
		ASSERT((i>=m_iMin)&&(i<=m_iMax));
		return m_array + i;
	};

	inline const type* ptr() const
	{
		return m_array;
	};

	inline type* ptr()
	{
		return m_array;
	};

	int length() const{return m_length;};
protected:
	int m_length;
	type* m_array;
	int m_iMin;
	int m_iMax;
	bool m_OwnsMemory;
};

template<class T> class CVarArray : public Array<T,T&>
{
public:
       CVarArray(int theIncrStep = 256):Array<T,T&>(0),
	m_size(0),m_incrStep(theIncrStep)
	{
		ASSERT(theIncrStep > 0);
	}

    /*
    CVarArray(const Array<T,T&>& array, int theIncrStep):Array<T,T&>(0),
    m_size(0),m_incrStep(theIncrStep)
       {
       	ASSERT(theIncrStep > 0);

       	int size = array.length()/theIncrStep;
       	size *= (theIncrStep+1);

		Array<T,T&>::m_array = (T*)Calloc(size, sizeof(T));//allocating bigger array
		if (!Array<T,T&>::m_array) throw "not enough memory";
		m_size = size;
		Array<T,T&>::m_iMin = 0;
		Array<T,T&>::m_iMax = array.length() - 1;
		Array<T,T&>::m_length = array.length();
       	setIndexShift(array.getIndexShift());
       	Array<T,T&>::copy(array);
       }*/

	virtual void create(int size)
	{
		int shift = -Array<T,T&>::m_iMin;
		Array<T,T&>::m_iMin = 0;
		Array<T,T&>::m_iMax = -1;

		Array<T,T&>::create(size);
		m_size = size;
		Array<T,T&>::setIndexShift(shift);
	}

	void moveIndex(int theShift)
	{
		if(m_size>0){
			Array<T,T&>::m_array+=theShift;
		}
		Array<T,T&>::m_iMin-=theShift;
		Array<T,T&>::m_iMax-=theShift;
	};
		
	void setIndexShift(int theShift)
	{
		if(m_size>0)
		{
			// moving to shift 0
			Array<T,T&>::m_array += Array<T,T&>::m_iMin;
		}
				
		// adjusting shift
		Array<T,T&>::m_iMax = Array<T,T&>::m_length - 1 - theShift;
		Array<T,T&>::m_iMin = -theShift;

		if(m_size>0)
		{
			// moving to new shift
			Array<T,T&>::m_array -= Array<T,T&>::m_iMin;
		}
	};

	int add(const T& theNewElement)
	{
		if (m_size<=0) {
			//create(m_incrStep);
			Array<T,T&>::m_array = (T*)Calloc(m_incrStep, sizeof(T));//allocating bigger array
			if (!Array<T,T&>::m_array) return 0;
			m_size = m_incrStep;

			Array<T,T&>::m_iMax = Array<T,T&>::m_length - 1 + Array<T,T&>::m_iMin;
			Array<T,T&>::m_array -= Array<T,T&>::m_iMin;
				
			return add(theNewElement);
		};
		if(Array<T,T&>::m_length==m_size)//memory for the array must be reallocated
		{
			T* newArray = (T*)realloc(Array<T,T&>::m_array + Array<T,T&>::m_iMin,(m_size+m_incrStep)*sizeof(T));//allocating bigger array
			if(!newArray)
				return 0;
			Array<T,T&>::m_array = newArray-Array<T,T&>::m_iMin;
			m_size += m_incrStep;
		}
		Array<T,T&>::m_length++;
		Array<T,T&>::m_iMax = Array<T,T&>::m_length - 1 + Array<T,T&>::m_iMin;
		Array<T,T&>::m_array[Array<T,T&>::m_iMax] = theNewElement;
		return 1;
	}

	T remove(int aIndex){
		T result = (*this)[aIndex];//ASSERT is done inside this statement
		if (aIndex!=Array<T,T&>::m_iMax)
			memmove(Array<T,T&>::m_array+aIndex, Array<T,T&>::m_array+aIndex+1, (Array<T,T&>::m_iMax-aIndex)*sizeof(T));
		Array<T,T&>::m_iMax--;
		Array<T,T&>::m_length--;
		return result;
	}

	virtual void reset()
	{
		if (Array<T,T&>::m_length) {
		Array<T,T&>::setIndexShift(0);
		Array<T,T&>::m_length = 0;
		}
	}
protected:
	int m_size;
	int m_incrStep;
};

typedef CVarArray<double> CVarVector;

template <class pointer, class pointee>
class CBasePointerArray :public CVarArray<pointer>
{
public:
	CBasePointerArray(int size=0)
	{
		if(size>0)
			CVarArray<pointer>::create(size);
	};
	pointee& operator()(int i)
	{
		ASSERT((*this)[i]!=NULL);
		return (*((*this)[i]));
	};
	const pointee& operator() (int i) const
	{
		ASSERT((*this)[i]!=NULL);
		return (*((*this)[i]));
	};
};

template <class type>
class CSmartPointerArray : public CBasePointerArray<SmartPtr<type>, type>{
public:
	CSmartPointerArray(int size=0) : CBasePointerArray<SmartPtr<type>, type>(size){};
	virtual void reset()
	{
		for(int i=CBasePointerArray<SmartPtr<type>, type>::m_iMin; i<=CBasePointerArray<SmartPtr<type>, type>::m_iMax; i++)
			(*this)[i] = NULL;
		CBasePointerArray<SmartPtr<type>, type>::reset();
	}
};

template <class type>
class CPointerArray :public CBasePointerArray<type*, type>
{
public:
	CPointerArray(int size=0) : CBasePointerArray<type*, type>(size){};

	void Delete(int aIndex){
		delete (CVarArray<type*>::remove(aIndex));
	}

	void DeleteAll(){
		for(int i=CPointerArray<type>::m_iMin;i<=CPointerArray<type>::m_iMax;i++){
			delete (*this)[i];
			(*this)[i] = NULL;
		}
		CPointerArray<type>::m_length = 0;
	}

};

template<class type>
class CAutoDeletePtrArray  : public CPointerArray<type>
{
public:
	CAutoDeletePtrArray(int size=0):CPointerArray<type>(size){};
	virtual ~CAutoDeletePtrArray()
	{
		for(int i=CPointerArray<type>::m_iMin;i<=CPointerArray<type>::m_iMax;i++)
			delete (*this)[i];
	};
};

typedef Array<double,double&> CBaseVector;

class CVector : public CBaseVector
{
public:
	CVector(double* aArray, int aLength):CBaseVector(aArray, aLength){};
	CVector(int size=0):CBaseVector(size){};
	CVector(const CBaseVector& array):CBaseVector(array){};
	CVector(const std::vector<double>& std_vec):CBaseVector(std_vec){};
#ifdef USE_GSL
	CVector(const gsl_vector *v) { create_gsl(v); }
	gsl_vector* GslVector() const;
	void create_gsl(const gsl_vector *v);
	void copy_gsl(const gsl_vector *v);
#endif
	CVector& operator=(const CBaseVector& theVector);
	CVector& operator+=(const CBaseVector& theVector);
	CVector& operator-=(const CBaseVector& theVector);
	CVector& operator*=(double theMult);
	CVector& operator/=(double theMult);
	CVector& operator+=(const int i)
	{
		moveIndex(i);
		return *this;
	};

	CVector& operator++()
	{
		moveIndex(1);
		return *this;
	}

	CVector& operator-=(const int i)
	{
		moveIndex(-i);
		return *this;
	};

	CVector& operator--()
	{
		moveIndex(-1);
		return *this;
	}

	//set all vector components to zero
	void reset();

	//find the nearest to xValue x-array element (increasing x-array supposed)
	//if value is out of range, returns -1
	static int findX(double xValue,double* xArray,int xArraySize);
	//find the nearest to xValue x-array element which is less than the argument (increasing x-array supposed)
	//if value is out of range, returns -1
	static int findLeftX(double xValue,double* xArray,int xArraySize);

	//find the nearest to xValue x-array element (increasing x-array supposed)
	//if value is out of range, returns <minimal array index>-1
	int findX(double xValue) const;

	//find the nearest to xValue x-array element which is less than the argument (increasing x-array supposed)
	//if value is out of range, returns <minimal array index>-1
	int findLeftX(double xValue) const;
};

typedef CVector* PVector;

class CMatrix
{
public:
	enum EMatrixType
	{
		rectangleE = 0,
		squareE,
		triangleLeftE,
		triangleRightE
	};
	CMatrix();
	//read rectangle matrix from a stream
	CMatrix(std::istream& aInput);
	CMatrix(int iSize, int jSize);
	CMatrix(int size, EMatrixType type = squareE); 
	void create(int iSize, int jSize);
	void create(int size, EMatrixType type = squareE);
	void reset();
	//This function solves the system A x = b in-place using Householder transformations (gsl_linalg_HH_svx)
	//where A - this matrix, b-the function argument
	//On input x should contain the right-hand side b, which is replaced by the solution on output.
	void solveLinearSystem (CVector& b) const;

	virtual ~CMatrix();
	inline CVector& operator[](int i){return *(m_array[i]);};
	inline const CVector& operator[](int i)const{return *(m_array[i]);};
	inline EMatrixType getMatrixType(){return m_type;};
	inline void setFirstIndexShift(int i){m_array.setIndexShift(i);};
	inline int SizeI() const {return m_array.length();}
	void multiply(double aValue);
protected:
	Array<PVector,PVector&> m_array;
	EMatrixType m_type;
};

class CMatrixAddOnlyView
{
public:
	inline CMatrixAddOnlyView(CMatrix& aMatrix) : fMatrix(aMatrix){}
	inline void Add(int i, int j, double aValue) {fMatrix[i][j] += aValue;}
private:
	CMatrix& fMatrix;
};

typedef CMatrix* PMatrix;


template<class type> // operations > and < must be defined for type
class COrderedList
{
public:
	COrderedList(int maxSize):
	m_array(maxSize),
	m_size(0)
	{
		m_order = new int[maxSize];
	};
	~COrderedList()
	{
		delete[] m_order;
	};
		
	int add(type element)
	{
		ASSERT(m_size<m_array.length());
		int no=searchOrder(element);
		if(no!=m_size)
			memmove(m_order+no+1, m_order+no, (m_size-no)*sizeof(int));
		m_order[no] = m_size;
		m_array[m_size] = element;
		m_size++;
		return no;
	};
	
	type& operator[](int i)
	{
		ASSERT(i>=0 && i<m_size);
		return m_array[m_order[i]];
	};

	//returns i-th element original order (order in which it was added to the list)
	int getOriginalOrder(int i)
	{
		return m_order[i];
	};
	int length(){return m_size;}
	double getMinValue(){return (*this)[0];}
	double getMaxValue(){return (*this)[m_size-1];}
	
protected:

	// returns number which will be assigned to the element being added
	int searchOrder(type& element)
	{
		if((!m_size) || element < m_array[m_order[0]])
			return 0;
			
		int iLeft = 0;
		int iRight = m_size;
		
		while(iRight-iLeft>1)
		{	
			int i=(iLeft + iRight)/2;
			if(element > m_array[m_order[i]])
				iLeft = i;
			else
				iRight = i;
		};
		return iRight;
	};
	Array<type, type&>	m_array;
	int*			m_order;
	int			m_size;
};

#endif // !defined(AFX_VECTOR_H__EA1E0B20_B7EC_11D5_885E_444553540001__INCLUDED_)
