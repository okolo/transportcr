// TableFunc.h
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TABLEFUNC_H__9C824721_C230_11D5_885E_444553540001__INCLUDED_)
#define AFX_TABLEFUNC_H__9C824721_C230_11D5_885E_444553540001__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include "Ranges.h"
#include <math.h>
#include "FilePtr.h"
#include "Function.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


class CTableReader;


class CTableFunction : public IFunction
{
public:
	CTableFunction(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.);
	CTableFunction(const CRanges& _ranges, CVector& _y);//takes x from CRanges::m_midEnergies

	//takes ownership of _table
	CTableFunction(CTableReader* _table, double _leftVal=0., double _rightVal=0.);
	
	virtual ~CTableFunction();

	//integrate f(x)g(x) from x_min to x_max 
	virtual double integrate(double g(double)=unitFunc) const;
	virtual double integrate(const IFunction& aFunc) const;

	/*
		calculates approximate integral this(t)g(t)dt from m_x[0] to x, m_x[0] < x < m_x[m_x.length-1]
		not implemented in base class
	*/
	virtual IFunction* integral(double g(double)=unitFunc) const;
	virtual IFunction* integral(const IFunction& aFunc) const;

	void init();
	virtual void print(const char* fileName) const;

	/// also extrapolates values of function beyond table range
	virtual double fApprox(double _x) const = 0;

	//get derivative (calculate if needed)
	// use every <period> point for calculation
	// period>1 may be usefull for stochastic or nonaccurate functions
	virtual CTableFunction* getDerivative(int period=1) const;

	/*
		create instance of itself 
		The method is used by IFunction* CTableFunction::integral(double g(double)=unitFunc)
		method to create object being returned as result.
	*/
	virtual CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const = 0;

	void takeDataOwnership();

	inline double firstX() const { return m_x[0];}
	inline double lastX() const { return m_x[m_x.length()-1];}
	inline double firstF() const { return m_y[0];}
	inline double lastF() const { return m_y[m_y.length()-1];}
	inline void setLeftVal(double aLeftVal){ m_leftVal = aLeftVal;}
	inline void setRightVal(double aRightVal){ m_rightVal = aRightVal;}
protected:
	virtual void prepare()=0;
	const CVector& m_x;
	const CVector& m_y;
	double m_leftVal;
	double m_rightVal;
	bool m_isReady;
	bool  m_ownsVectors;
	CTableReader* m_table;
};

template <class type>
void testTableFunction(double func(double), double xMin, double xMax, bool isLogScale = false,
					   UINT accuracy = 100, const char* fileName=NULL/* PLT_DIR/Func.test*/)
{
	if (fileName==NULL)
	{
		fileName = "Func.test";
	}

	ASSERT(xMin < xMax);
	if(isLogScale)
	{
		ASSERT(xMin>0.);
	}

	CVector x(accuracy+1);
	CVector y(accuracy+1);
	double dx = isLogScale?(pow(xMax/xMin,1./accuracy)):((xMax - xMin)/accuracy);
	double xCurrent = xMin;
	UINT i=0;
	for(; i <= accuracy; i++)
	{
		x[i] = xCurrent;
		y[i] = func(xCurrent);
		if (isLogScale)
		{
			xCurrent*=dx;
		}else
		{
			xCurrent+=dx;
		}
	}
	type f(x,y);
	f.init();

	CFilePtr file(Fopen(fileName,"wt"));

	if (isLogScale)
	{
		dx = pow(dx, 1./10.);
	}
	else
	{
		dx /= 10.;
	}
	
	fprintf(file,"# x\tF_table\tF_real\n");
	for(xCurrent = xMin; xCurrent < xMax; )
	{
		fprintf(file, "%lg\t%lg\t%lg\n", xCurrent, f(xCurrent),func(xCurrent));
		if (isLogScale)
		{
			xCurrent*=dx;
		}else
		{
			xCurrent+=dx;
		}
	}		
}

namespace PropagDepreciated
{
/**
Implements linear approximation of table function in logarifmic scale
*/
class CLogScaleFunc  : public CTableFunction
{
public:
	CLogScaleFunc(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.):CTableFunction(_x, _y, _leftVal, _rightVal){init();};
	CLogScaleFunc(const CRanges& _ranges, CVector& _y):CTableFunction(_ranges, _y){init();};//takes x from CRanges::m_midEnergies
	CLogScaleFunc(CTableReader* _table, double _leftVal=0., double _rightVal=0.):CTableFunction(_table, _leftVal, _rightVal){init();};
	
	virtual ~CLogScaleFunc();
	double f(double _x) const;
	double fApprox(double _x) const;

	CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const{return new CLogScaleFunc(_x, _y, _leftVal, _rightVal);};

protected:
	virtual void prepare();
	CVector m_alphas;

};
}// end PropagDepreciated


/**
Cubic spline interpolation
*/
class CCSplineFunc  : public CTableFunction
{
public:
	CCSplineFunc(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.):CTableFunction(_x, _y, _leftVal, _rightVal){init();};
	CCSplineFunc(const CRanges& _ranges, CVector& _y):CTableFunction(_ranges, _y){init();};//takes x from CRanges::m_midEnergies
	CCSplineFunc(CTableReader* _table, double _leftVal=0., double _rightVal=0.):CTableFunction(_table, _leftVal, _rightVal){init();};
	
	virtual ~CCSplineFunc();
	inline void setXminSpline(double aXmin){ m_Xmin=aXmin;};
	inline void setXmaxSpline(double aXmax){ m_Xmax=aXmax;};



	/// CTableFunction overridden methods
	virtual double f(double _x) const;

	virtual double fApprox(double _x) const;

	//get derivative (calculate if needed)
	// use every <period> point for calculation
	// period parameter is not used in CCSplineFunc::getDerivative implementation
	virtual CTableFunction* getDerivative(int period=1) const;

	/*
		get derivative (calculate if needed)
	*/
	CTableFunction* get2ndDerivative() const;

	/*
		create instance of itself 
		The method is used by IFunction* CTableFunction::integral(double g(double)=unitFunc)
		method to create object being returned as result.
	*/
	CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const{return new CCSplineFunc(_x, _y, _leftVal, _rightVal);};


protected:
	virtual void prepare();
	gsl_interp_accel*		m_accel;
	gsl_spline*				m_spline;

	/// spline validity range (default values: m_x[0] and m_x[m_x.length()-1])
	double					m_Xmin;
	double					m_Xmax;
};

///// Uses cubic spline of log(F), note that zero and negative points are just ignored,
///// error messages are not produced in release configuration!!!
///// to varify that input data is suitable for log cspline fit call assertData()
class CLogCSplineFunc  : public CTableFunction
{
public:
	CLogCSplineFunc(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.):CTableFunction(_x, _y, _leftVal, _rightVal){init();};
	CLogCSplineFunc(const CRanges& _ranges, CVector& _y):CTableFunction(_ranges, _y){init();};//takes x from CRanges::m_midEnergies
	CLogCSplineFunc(CTableReader* _table, double _leftVal=0., double _rightVal=0.):CTableFunction(_table, _leftVal, _rightVal){init();};

	/// CTableFunction overridden methods
	virtual double f(double _x) const;

	virtual double fApprox(double _x) const;

	void assertData();
	CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const{return new CLogCSplineFunc(_x, _y, _leftVal, _rightVal);};

protected:
	virtual void prepare();
	gsl_interp_accel*		m_accel;
	gsl_spline*				m_spline;

	CVarVector				m_logF;
	CVarVector				m_adjustedX;

	/// spline validity range (default values: first and last X with nonzero F)
	double					m_Xmin;
	double					m_Xmax;
};

class MatrixFunction : public IFunction2
{
public:
	MatrixFunction(const std::string& aFile);
	double f(double x, double y) const;
	virtual double MinArg(int aArgNo) const {return aArgNo==1?xMin:yMin;};
	virtual double MaxArg(int aArgNo) const {return aArgNo==1?xMax:yMax;};
private:
	SafePtr<CMatrix> fData;
	double xMin;
	double yMin;
	double xMax;
	double yMax;
	double xStep;
	double yStep;
	bool logScaleX;
	bool logScaleY;
};

enum EExtensionRule //below xmin and/or above xmax
{
	ExtNone = 0,//cut to zero
	ExtConst = 1,//continue as constant
	ExtLinear = 2//continue as linear function (in log or linear scale)
};

/**
Implements simple lenear approximation of table function
*/
class CLinearFunc  : public CTableFunction
{
public:
	CLinearFunc(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.):
	  CTableFunction(_x, _y, _leftVal, _rightVal),
	  m_leftDerivative(0.),
	  m_rightDerivative(0.)
	  {
		  init();
	  };
	CLinearFunc(const CRanges& _ranges, CVector& _y):
	  CTableFunction(_ranges, _y),
	  m_leftDerivative(0.),
	  m_rightDerivative(0.)
	  {
		  init();
	  };//takes x from CRanges::m_midEnergies

	//takes ownership of _table
	CLinearFunc(CTableReader* _table, double _leftVal=0., double _rightVal=0.):
	  CTableFunction(_table, _leftVal, _rightVal),
	  m_leftDerivative(0.),
	  m_rightDerivative(0.)
	  {
		  init();
	  };

    void SetExtension(EExtensionRule aLimits, bool aLeft = true, bool aRight = true);

	double f(double _x) const;
	double fApprox(double _x) const;

	/*
		create instance of itself 
		The method is used by IFunction* CTableFunction::integral(double g(double)=unitFunc)
		method to create object being returned as result.
	*/
	CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const{return new CLinearFunc(_x, _y, _leftVal, _rightVal);};
	
protected:
	virtual void prepare();
// fields 
//	CVector m_derivativeX;
//	CVector m_derivativeY;
private:
	double m_leftDerivative;
	double m_rightDerivative;
};

class CLogScaleLinearFunc  : public CTableFunction
{
public:
	CLogScaleLinearFunc(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.):
	  CTableFunction(_x, _y, _leftVal, _rightVal)
	  {
		  init();
	  };
	CLogScaleLinearFunc(const CRanges& _ranges, CVector& _y):
	  CTableFunction(_ranges, _y)
	  {
		  init();
	  };//takes x from CRanges::m_midEnergies

	//takes ownership of _table
	CLogScaleLinearFunc(CTableReader* _table, double _leftVal=0., double _rightVal=0.):
	  CTableFunction(_table, _leftVal, _rightVal)
	  {
		  init();
	  };

    void SetExtension(EExtensionRule aLimits, bool aLeft = true, bool aRight = true);

	double f(double _x) const;
	double fApprox(double _x) const;

	/*
		create instance of itself
		The method is used by IFunction* CTableFunction::integral(double g(double)=unitFunc)
		method to create object being returned as result.
	*/
	CTableFunction* createInstance(CVector& _x, CVector& _y, double _leftVal=0., double _rightVal=0.) const{return new CLogScaleLinearFunc(_x, _y, _leftVal, _rightVal);};

protected:
	virtual void prepare();
private:
	SafePtr<CLinearFunc> fLogXY;
	static const double Log0;
};

//typedef CCSplineFunc CDefaultTableFunc;
typedef CLinearFunc CDefaultTableFunc;
//typedef PropagDepreciated::CLogScaleFunc CDefaultTableFunc;

#endif // !defined(AFX_TABLEFUNC_H__9C824721_C230_11D5_885E_444553540001__INCLUDED_)
