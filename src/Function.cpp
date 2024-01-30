
#include <math.h>
#include "Function.h"

double unitFunc(double)
{
	return 1.;
}

double unitFunc2(double,double)
{
	return 1.;
}

double xFunc(double x)
{
	return x;
}

CFunc gUnitFunc(unitFunc);
CFunc gLinearFunc(xFunc);
CFunc2 gUnitFunc2(unitFunc2);

double CFunc::integrate(double xMin, double xMax, double g(double), int ac)
{
	CFunc func(g);
	return FuncUtils::integrate(xMin, xMax, this, &func, ac);
}

double CFunc::integrate(double xMin, double xMax, const IFunction* g, int ac)
{
	return FuncUtils::integrate(xMin, xMax, this, g, ac);
}

double FuncUtils::integrate(double xMin, double xMax, const IFunction* g1, const IFunction* g2, int ac)
{
	double dx=(xMax-xMin)/ac;
	int i=0;
	double result = 0.;
	for (double x=xMin+0.5*dx; i<ac; x+=dx,i++) {
		result += (*g1)(x)*(*g2)(x);
	}
	result *= dx;
	return result;
}

double FuncUtils::integrateLogscale(double xMin, double xMax, const IFunction* g, int ac)
{
	ASSERT(xMin > 0);
	ASSERT(xMax > 0);
	double dx = log(xMax/xMin)/ac;
	double mult = exp(dx);
	int i=0;
	double result = 0.;
	for (double x=xMin*sqrt(mult); i<ac; x*=mult,i++) {
		result += (((*g)(x)) * x);
	}
	result *= dx;
	return result;
}

double CPowerLawFunction::f(double _x) const
{
	ASSERT(_x>=0.);
	return pow(_x, m_power);
}

double CPowerLawFunction::f_inverse(double _x) const
{
	ASSERT(_x>=0.);
	ASSERT(m_power!=0.);
	return pow(_x, 1./m_power);
}

IFunction* CPowerLawFunction::integral(double xMin, double xMax, double g(double))
{
	if (g==unitFunc) {
		return integralPower(xMin, 0.);
	}
	else if(g==xFunc){
		return integralPower(xMin, 1.);
	}

	NOT_IMPLEMENTED
	return NULL;
}

/*
Implementation of power law function integration
*/
class CPowerLawFunctionIntegral : public IFunction{
public:
	CPowerLawFunctionIntegral(double xMin, double power):
	m_add(0.),
	m_power(power+1.)
	{
		ASSERT(m_power!=0.);
		m_mult = 1./(m_power);
		m_add = -f(xMin);
	}
	double f(double _x) const
	{
		return m_mult*pow(_x, m_power) + m_add;
	}
	double f_inverse(double _x) const
	{
		return pow((_x-m_add)/m_mult,1./m_power);
	}
protected:
	double m_mult;
	double m_add;
	double m_power;
};
class CPowerLawLogFunctionIntegral : public IFunction{
public:
	CPowerLawLogFunctionIntegral(double xMin):
	m_xMin(xMin)
	{
		ASSERT(m_xMin>0.);	
	}
	double f(double _x) const
	{
		return log(_x/m_xMin);
	}
	double f_inverse(double _x) const
	{
		return m_xMin*exp(_x);
	}
protected:
	double m_xMin;
};

IFunction* CPowerLawFunction::integralPower(double xMin, double gPower)
{
	double power = m_power+gPower;
	if (power==-1.) {
		return new CPowerLawLogFunctionIntegral(xMin);
	}
	return new CPowerLawFunctionIntegral(xMin, m_power+gPower);
}

static double GslProxyFunc (double x, void * params)
{
	const IFunction* f = (const IFunction*)params;
	return f->f(x);
}

FuncUtils funcUtils;

FuncUtils::~FuncUtils()
{
	if(gslQAGintegrator)
		gsl_integration_workspace_free (gslQAGintegrator);
}

double FuncUtils::gsl_integration_qag (const IFunction* aFunc, double aXmin, double aXmax, double epsabs, double epsrel, size_t limit, int key)
{
	if(gslQAGintegrator==0)
		gslQAGintegrator = gsl_integration_workspace_alloc (limit);
	else if(gslQAGintegrator->limit < limit)
	{
		gsl_integration_workspace_free (gslQAGintegrator);
		gslQAGintegrator = gsl_integration_workspace_alloc (limit);
	}

	gsl_function f;
	f.params = (void*)aFunc;
	f.function = GslProxyFunc;
	double result, abserr;
	int failed = ::gsl_integration_qag (&f, aXmin, aXmax, epsabs, epsrel, limit, key, gslQAGintegrator, &result, &abserr);
	if(failed)
	{
		ASSERT(0);
		ThrowError("Integration failed with code " + ToString(failed));
	}
	return result;
}
