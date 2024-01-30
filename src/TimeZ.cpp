// TimeZ.cpp: implementation of the CTimeZ class.
//
//////////////////////////////////////////////////////////////////////

#include "TimeZ.h"
#include <math.h>
#include "Addfunc.h"
#include "Vector.h"
#include "TableFunc.h"
#include "Ranges.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include "InjectionSpectra.h"
#include "Units.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double CTimeZ::Lv = 0.;
double CTimeZ::Lm = 0.;
double CTimeZ::t0 = 0.;
double CTimeZ::sqrtLm = 0.;
double CTimeZ::sqrtLv = 0.;
double CTimeZ::Hubble = 0.;

const int CTimeZ::Accuracy = 100000;

static CTimeZ s_timeZ;

void CTimeZ::setLv(double _Lv)
{
  //ASSERT((_Lv>=0.)&&(_Lv<=1.));
	if(_Lv<0)
	  _Lv = 0;
	Lv = _Lv;
	Lm = 1.-Lv;	
	sqrtLv = sqrt(Lv);
	sqrtLm = sqrt(Lm);
}

void CTimeZ::init(double _Lv, double _Hubble_in_km_s_Mpc)
{
	Hubble = _Hubble_in_km_s_Mpc;
	Hubble /= (units.Mpc_cm *1e-5);/* Hubble constant in 1/s  */
	Hubble *= units.Tunit;

	if((_Lv>=0)&&(Lv<=1.))
		setLv(_Lv);
	if(Lv <= 0)//Lyambda = 0
	{
		t0 = 2./3./Hubble;
	}
	else
	{
		t0 = log((2. + 2.*sqrtLv - Lm)/Lm)/(3.*sqrtLv)/Hubble;
		//t0 = 2./3./H/sqrtLv*log((1+sqrtLv)/sqrtLm);
	}
}

double CTimeZ::t2z(double _t)
 {
	 /*
 	if(Lv <= 0)//Lyambda = 0
 	{
		return pow(1.5*H*_t,-2./3.)-1.;
	}
	double c = exp(1.5*H*sqrtLv*_t)*sqrtLm;
	double x = 0.5*(c - Lm/c);
	return pow(Lv/x/x,1./3.)-1.;*/
 	if(!s_timeZ.t2zFunc)
		s_timeZ.initZ2T();
	return s_timeZ.t2zFunc->f(t0-_t);
}

double CTimeZ::z2t(double _z)
{
	/*
	if(Lv <= 0)//Lyambda = 0
	{
		return 2./3./H/pow(1.+_z,1.5);
	}

	double x = sqrtLv/pow(1.+_z,1.5);
	return 2./3./H/sqrtLv*log((x + sqrt(Lm + x*x))/sqrtLm);*/

	if(!s_timeZ.z2tFunc)
		s_timeZ.initZ2T();
	return t0 - s_timeZ.z2tFunc->f(_z);
}

double CTimeZ::diffD(double a)
{
	return 1./( Hubble*sqrt(a*(Lm+Lv*a*a*a)) );
}

//dDc/dz ( Dc is comoving distance as defined in astro-ph/9905116 )
double CTimeZ::diffDc(double aZ)
{
	return 1./H(aZ);
}

//Hubble parameter at epoch z
double CTimeZ::H(double _z)
{
	double a = 1.+_z;
	return Hubble*sqrt(Lm+Lv*a*a*a);
}

void CTimeZ::initZ2D()
{
	ASSERT(Ranges().Zmax()>0.);// must be initialized
	ASSERT(BC().ss1>1.); // must be initialized
	double zMax = Ranges().Zmax() + 100.;
	
	z = new CVector(Accuracy+1);
	d = new CVector(Accuracy+1);
	(*z)[0]=0;
	(*d)[0]=0;
	
	double aMin = 1./(zMax+1)/BC().ss1;
	double da = (1.-aMin)/Accuracy;
	double da2 = 0.5*da;
	double a = 1.-da;
	double prev = diffD(1.);
	double sum = 0.;

	for(int i=1; i<=Accuracy; i++, a -= da)
	{
		sum += prev;
		sum += 4.*diffD(a+da2);
		prev = diffD(a);
		sum += prev;
		(*z)[i] = 1./a-1.;
		(*d)[i] = da*sum/6.;
	}
	z2dFunc = new CLinearFunc(*z,*d,0,(*d)[Accuracy]);
	d2zFunc = new CLinearFunc(*d,*z,0,(*z)[Accuracy]);
}

double CTimeZ::z2d(double _z)
{
	if(!s_timeZ.z2dFunc)
		s_timeZ.initZ2D();
	return s_timeZ.z2dFunc->f(_z);
}

double CTimeZ::d2z(double _d)
{
	if(!s_timeZ.d2zFunc)
		s_timeZ.initZ2D();
	return s_timeZ.d2zFunc->f(_d);
}

void CTimeZ::initZ2Dc()
{
	ASSERT(Ranges().Zmax()>0.);// must be initialized
	ASSERT(BC().ss1>1.); // must be initialized
	double zMax = Ranges().Zmax() + 100.;

	zDc = new CVector(Accuracy+1);
	dDc = new CVector(Accuracy+1);
	(*zDc)[0]=0;
	(*dDc)[0]=0;

	double dz = zMax/Accuracy;
	double dz2 = 0.5*dz;
	double z = dz;
	double prev = diffDc(0.);
	double sum = 0.;

	for(int i=1; i<=Accuracy; i++, z += dz)
	{
		sum += prev;
		sum += 4.*diffDc(z-dz2);
		prev = diffDc(z);
		sum += prev;
		(*zDc)[i] = z;
		(*dDc)[i] = dz*sum/6.;
	}

	z2DcFunc = new CLinearFunc(*zDc,*dDc,0,(*dDc)[Accuracy]);
}

double CTimeZ::z2Dc(double _z)
{
	if(!s_timeZ.z2DcFunc)
		s_timeZ.initZ2Dc();
	return s_timeZ.z2DcFunc->f(_z);
}

double CTimeZ::diffT(double aZ)
{
	double z1 = 1.+aZ;
	return 1./Hubble/z1/sqrt(Lm*z1*z1*z1+Lv);
}

void CTimeZ::initZ2T()
{
	ASSERT(Ranges().Zmax()>0.);// must be initialized
	ASSERT(BC().ss1>1.); // must be initialized
	double zMax = Ranges().Zmax() + 100.;
	zt = new CVector(Accuracy+1);
	tt = new CVector(Accuracy+1);
	(*zt)[0]=0;
	(*zt)[0]=0;
	
	double dz = zMax/Accuracy;
	double dz2 = 0.5*dz;
	double z = dz;
	double prev = diffT(0.);
	double sum = 0.;

	for(int i=1; i<=Accuracy; i++, z += dz)
	{
		sum += prev;
		sum += 4.*diffT(z-dz2);
		prev = diffT(z);
		sum += prev;
		(*zt)[i] = z;
		(*tt)[i] = dz*sum/6.;
	}
	z2tFunc = new CLinearFunc(*zt,*tt,0,(*tt)[Accuracy]);
	t2zFunc = new CLinearFunc(*tt,*zt,0,(*zt)[Accuracy]);
}

///The function is used to test redshift evolution in special case of the source with no evolution in comoving frame and alpha = 2
double CTimeZ::TestIntegral(double r /*r===min(Emax/E, 1+Zmax)*/)
{
	return TestIntegralF(r)-TestIntegralF(1.);
}

double Hy2F1(double a, double b, double c, double z)
{
	if(fabs(z)<=1)
		return gsl_sf_hyperg_2F1(a,b,c,z);
	z = 1./z;
	double mult1 = gsl_sf_gamma(b-a)*gsl_sf_gamma(c)/gsl_sf_gamma(b)/gsl_sf_gamma(c-a);
	mult1 *= pow(-z, a);
	double mult2 = gsl_sf_gamma(a-b)*gsl_sf_gamma(c)/gsl_sf_gamma(a)/gsl_sf_gamma(c-b);
	mult2 *= pow(-z, b);
	double f1 = gsl_sf_hyperg_2F1(a,a-c+1,a-b+1,z);
	double f2 = gsl_sf_hyperg_2F1(b,b-c+1,b-a+1,z);
	double result = mult1*f1 + mult2*f2;
	return result;
}

///The function is used to test redshift evolution in special case of the source with no evolution in comoving frame and alpha = 2
double CTimeZ::TestIntegralF(double a/*a===1+z*/)//F[1+z, alpha] = (1+z)^(-alpha)*(Lm*(1+z)^3 + Lv)^(-1/2)
{//TestIntegralF === Integrate[F[1+z, 1], z]
	double a3 = a*a*a;
	double arg = -(a3*Lm)/Lv;
	double Hypergeometric2F1 = Hy2F1(0.6666666666666666,0.5,1.6666666666666667,
         arg);
	double result = (sqrt(1 + (-1 + a3)*Lm)*
     (4*(-1 + Lm)*sqrt((-1 + Lm - a3*Lm)/
          (-1 + Lm)) + 
       a3*Lm*Hypergeometric2F1))/(4.*a*Lv*Lv*
     sqrt((-1 + Lm - a3*Lm)/(-1 + Lm)));
	return result;
}

void CTimeZ::SuppressionTest()
{
	double zCalculationMax = Ranges().Zmax();
	cerr << "redshift spectrum suppression test for m-alpha = -2 and Zmax = " << zCalculationMax << ":\n";
//here we simulate how spectrum should be suppressed at highest energies due to redshifting
//The analytic estimate was made for spectrum Q(E,z) = E^-alpha * (1+z)^(3+m) * Theta(Emax-E), for m-alpha = -2
	double norm = 1./TestIntegral(1.+zCalculationMax);
	for(double r=1; r<zCalculationMax + 1; r+=0.05)
		cerr << r << "\t" << TestIntegral(r)*norm << "\n";
}

Redshift::Redshift():
m_z(-1.)
{

}

Redshift::Redshift(double aZ)
{
	setZ(aZ);
}

void Redshift::setZ(double aZ)
{
	ASSERT(aZ>=0.); 
	m_z = aZ;
	m_dl = (1.+m_z);
	m_dv = m_dl*m_dl*m_dl;
	m_t = CTimeZ::z2t(aZ);
}

Redshift redshift;
