//#include <math.h>
#ifndef FUNCTION_H_TIPA_ALREADY_INCLUDED
#define FUNCTION_H_TIPA_ALREADY_INCLUDED

#include "Addfunc.h"
#include <gsl/gsl_integration.h>
#include <limits>

extern double unitFunc(double); // f(x)=1
extern double xFunc(double); //  f(x)=x
extern class CFunc gUnitFunc; // f(x)=1
extern class CFunc gLinearFunc;//  f(x)=x
extern class CFunc2 gUnitFunc2;//f(x,y)=1

class IFunction
{
public:
	virtual double f(double _x) const =0;
	inline double operator()(double _x) const {return f(_x);};
	virtual ~IFunction(){};

	/*
		returns F(x)=Integrate{f(t)g(t),{t,xMin,x}} for xMin < x < xMax
		the result must be deleted manually
	*/
	virtual IFunction* integral(double xMin, double xMax, double g(double)=unitFunc){NOT_IMPLEMENTED;return NULL;};

	// inverse function (if exists)
	virtual double f_inverse(double _x) const {NOT_IMPLEMENTED;return 0.;};
};

class IFunction2
{
public:
	virtual double f(double x, double y) const = 0;
	virtual double MinArg(int aArgNo) const {return std::numeric_limits<double>::min();}
	virtual double MaxArg(int aArgNo) const  {return std::numeric_limits<double>::max();}
	inline double operator()(double _x, double _y) const {return f(_x,_y);};
	virtual ~IFunction2(){};
};

class CFunc : public IFunction
{
public:
	CFunc(double _f(double)):m_f(_f){}
	double f(double _x) const{return m_f(_x);}
	double integrate(double xMin, double xMax, double g(double)=unitFunc, int ac=100);
	double integrate(double xMin, double xMax, const IFunction* g = &gUnitFunc, int ac=100);

protected:
	double (*m_f)(double);
};

class CFunc2 : public IFunction2
{
public:
	CFunc2(double _f(double, double)):m_f(_f){}
	double f(double _x, double _y) const{return m_f(_x,_y);}
protected:
	double (*m_f)(double, double);
};

class FuncUtils
{
public:
	FuncUtils():gslQAGintegrator(0){}
	~FuncUtils();
	static double integrateLogscale(double xMin, double xMax, const IFunction* g, int ac=100);
	static double integrate(double xMin, double xMax, const IFunction* g1, const IFunction* g2 = &gUnitFunc, int ac=100);

	inline static double integrateBode(double xMin, double xMax, const IFunction& g)
	{
		double leftVal = g(xMin);
		return (xMax-xMin)*integrateBodeFast(xMin, xMax, g, &leftVal);
	}

	//Bode's integration role, using left value and do not multiply by deltaX for faster integration
	inline static double integrateBodeFast(double xMin, double xMax, const IFunction& g, double* leftGvalue)
	{
		const double a1=0.07777777777777778;//  7/90;
		const double a2=0.3555555555555556;// 16/45;
		const double a3=0.1333333333333333;// 6/45;
		double x2=0.5*(1.5*xMin+0.5*xMax);
		double x3=0.5*(xMin + xMax);
		double x4=0.5*(0.5*xMin + 1.5*xMax);
		double result = *leftGvalue;
		*leftGvalue = g(xMax);
		result += *leftGvalue;
		result *= a1;
		result += (a2*(g(x2)+g(x4)) + a3*g(x3));
		return result;
	}

	//QAG adaptive integration
	double gsl_integration_qag (const IFunction* aFunc, double aXmin, double aXmax, double epsabs, double epsrel, size_t aMaxSubIntervals, int key=GSL_INTEG_GAUSS15);
private:
	gsl_integration_workspace* gslQAGintegrator;
};

extern FuncUtils funcUtils;

class TSuperPosition : public IFunction
{
public:
	TSuperPosition(const IFunction& aF1, const IFunction& aF2):iF1(aF1), iF2(aF2){};
	double f(double _x) const{return iF1(iF2(_x));}
private:
	const IFunction& iF1;
	const IFunction& iF2;
};

class TProduct : public IFunction
{
public:
	TProduct(const IFunction& aF1, const IFunction& aF2):iF1(aF1), iF2(aF2){};
	double f(double _x) const{return iF1(_x)*iF2(_x);}
private:
	const IFunction& iF1;
	const IFunction& iF2;
};

class CPowerLawFunction : public IFunction
{
public:
	CPowerLawFunction(double _power):m_power(_power){}
	double f(double _x) const;
	double f_inverse(double _x) const;
	IFunction* integral(double xMin, double xMax, double g(double)=unitFunc);

	/*
		returns F(x)=Integrate{f(t)*t^(gPower),{t,xMin,x}} for xMin < x < xMax
		the result must be deleted manually
	*/
	IFunction* integralPower(double xMin, double gPower);
protected:
	double m_power;
};

class ArgShiftFunc : public IFunction
{
public:
    //takes ownership of aFunc
    ArgShiftFunc(IFunction* aFunc, double aArgShift) : fF(aFunc), fDx(aArgShift)
            {

            }
    double f(double _x) const{
        return fF->f(_x + fDx);
    }

private:
    SafePtr<IFunction>  fF;
    double              fDx;
};

#endif //ifndef FUNCTION_H_TIPA_ALREADY_INCLUDED
