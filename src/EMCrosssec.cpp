/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888
888888888888            Cross section function definitions        88888888888888888888888
88888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/

#include <math.h>
#include "Addfunc.h"
#include "Ranges.h"
#include "gsl/gsl_sf_dilog.h"
#include "Units.h"
#include "EMCrosssec.h"

const double  EMCrosssec::C1=8.3647043703e-5;
        /* C1=0.5*Pi/137.035989561/137.035989561=
           =3/16*TompsCS*Me^2=Pi/2*alpha^2 */
const double EMCrosssec::ln2=0.693147180559945;

EMCrosssec::EMCrosssec(double aMass):
fMass(aMass),
fSquareMass(aMass*aMass)
{
}

#define ICS_JONES
#define PP_ZDZIARSKI
#define Li2_GSL

double EMCrosssec::li2(double aX)
{
	return  aX<1e-8 && aX>-1e-8 ?  aX + 0.25*aX*aX : gsl_sf_dilog (aX);
}

double EMCrosssec::intPICSph_Jones(double A, double q)
{
	double A2 = A*A;
	double A3 = A*A2;
	double A4 = A2*A2;
	double q2 = q*q;
	double q3 = q2*q;
	double Aq = A*q;
	double ln_Aq_p1 = (Aq>1e-5)?log(1. + Aq):(Aq*(1.+Aq*(-0.5+0.33333333333*Aq)));
	double Aq_p1 = 1. + Aq;
	double B = Aq*(2.*Aq+Aq);
	double Aq_p12_1 = B>1e-5 ? (1./(1.+B)) : (1.-B+B*B-B*B*B);

	double result = (-16*A*q - 10*A2*q + 2*A3*q - 
     24*A2*q2 - 13*A3*q2 + A4*q2 - 8*A3*q3 - 
     2*A4*q3 + 16*ln_Aq_p1 + 
     18*A*ln_Aq_p1 + 2*A2*ln_Aq_p1 + 
     32*A*q*ln_Aq_p1 + 
     36*A2*q*ln_Aq_p1 + 
     4*A3*q*ln_Aq_p1 + 
     16*A2*q2*ln_Aq_p1 + 
     18*A3*q2*ln_Aq_p1 + 
     2*A4*q2*ln_Aq_p1 + 
     8*A*Aq_p1*log(q)*
      (-Aq + Aq_p1*ln_Aq_p1))*Aq_p12_1/
   (4.*A3);

	result += (2.*li2(-Aq)/A2);

	ASSERT(result>=0);
	return result;
}

double EMCrosssec::intPICSph_Jones(double E, double b, double nE)//return intergal of P from 0 to nE
{
	double A = 4.*b*E/fSquareMass;
	double mult = C1 * 16.0 /fSquareMass;
	double r = nE/E;

	double q = 1.;
	if(r<1.)
	{
		if(r<0.001)
		{
			q = r*(1. + r*(1. + r*(1+r)))/A;
		}else
		{
			q = r/(1.-r)/A;
		}
	}

	if(q>1)
		q=1.;
	double result = intPICSph_Jones(A, q);
	ASSERT(result>=0.);
	return result*mult;
}

//maximal nE/E for secondary photon in ICS (see PICSphR(double Eb, double r))
double EMCrosssec::rMaxICSph(double Eb)
{
	double A = 4.*Eb/fSquareMass;
	return A/(1.+A);
}

/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  /*Inverse Compton Scattering angle averaged differential cross section
    (per photon energy) multiplied by E

	source: Formula (44) of Frank C. Jones, Phys Rev Volume 167, Number 5 (1968), pp. 1159-1166
	Eb = E*b
	r = nE/E;
	E-cosmic ray electron energy
	b-background photon energy
	nE-outgoing photon energy
	// the formula assumes that nE > b
	// this is lowest order term in double expansion in b/E and 1/gamma^2
*/
double EMCrosssec::PICSphR(double Eb, double r)
{
	double A = 4.*Eb/fSquareMass;
	double rMax = A/(1.+A);
	if(r>=rMax) return 0.;
	
	double eps = rMax-r;
	double k,q,q_1;//q_1 = 1-q
	if(A*eps<1e-3)
	{//calculating 1-q
		double A1eps=(1.+A)*eps;
		q_1=A1eps*(1-A1eps*(1+A1eps))*(1.+1./A);
		ASSERT(q_1<1 && q_1>0);
		k=(1.-q_1)*A;
		q = 1.-q_1;
	}
	else
	{
		if(r<0.001)
		{
			k = r*(1. + r*(1. + r*(1+r)));
		}else
		{
			k = r/(1.-r);
		}
		q = k/A;
		q_1 = 1.-q;
	}
	ASSERT(q<1);
	ASSERT(q>=0.);
	if(q==0.)
	{
		return C1 * 16.0 /(fSquareMass*A);	
	}
	double Log_q = (q_1<1e-5)? q*(-q_1*(q_1*(0.5 + q_1*(1./3. + 0.25*q_1)))) : log(q);
	double result = C1 * 16.0 /(fSquareMass*A) * (2.*q*Log_q + (1.+2.*q)*q_1 + 0.5 * k*k/(k+1.)*q_1);
	ASSERT_VALID_NO(result);
	return result;
}

double EMCrosssec::PICSph(double E,double b,double nE)
{
	double Eb = b*E;
	double r = nE/E;
	return PICSphR(Eb,r)/E;
}

double EMCrosssec::intPrdrICSph_Jones(double E, double b)
{//used for analythical estimate of mean secondary gamma energy
	double A = 4.*b*E/fSquareMass;
	double A2 = A*A;
	double A3 = A2*A;
	double A4 = A2*A2;
	double result = 0;
	if(A<1e-2)
	{
		result = 0.111111111111111 - 0.175*A + 0.245*A2 - 0.320634920634921*A3 + 0.399659863945578*A4;
	}
	else
	{
		double Ap12 = 1.+2.*A+A2;
		double lnA = log(A);
		double lnAp1 = log(1.+A);
		result = -(A*(72. + 156.*A + 96.*A2 + 11.*A3 + 4.*Ap12*Pi*Pi + 12.*Ap12*lnA*lnA)
			- 6.*Ap12*(12. + A*(12. + A))*lnAp1 + 24.*A*Ap12*li2(-1./A))/(12.*A4*Ap12);
	}
	ASSERT_VALID_NO(result);
	return result * C1 * 16.0 * A / fSquareMass;
}

double EMCrosssec::intPrdrICSph(double E, double b)
{//used for analythical estimate of mean secondary gamma energy from ICS
#ifdef ICS_JONES
	return intPrdrICSph_Jones(E,b);
#else
	ASSERT(false);
	throw "unsupported function call";
#endif
}

//source: Formula (47) of Frank C. Jones, Phys Rev Volume 167, Number 5 (1968), pp. 1159-1166
// see also F. C. Jones, Phys. Rev. 137, B1306 (1965), formulas (13) and (14)
//energy loss rate of electron 1/E dE/dt (should be integrated with gamma background)
double EMCrosssec::EnergyLossICS(double Eb)
{
	double a = 2.*Eb/fSquareMass;
	double a2 = a*a;
	double a3 = a2*a;

	double result = C1/fSquareMass*4./a2;
	double F = 0;
	if(a>1e-3)
	{
		double a21 = 1.+2.*a;
		F = (a+6.+3./a)*log(a21);
		F -= (22./3.*a3+24.*a2+18.*a+4.)/a21/a21;
		F -= 2.*(1.-li2(-2.*a));
	}
	else
	{//replacing by series
		F=8./9.*a*a*a-14./5.*a*a*a*a;
	}
	result *= F;
	ASSERT_VALID_NO(result);
	return result;
}

//energy loss rate of electron dE/dt (should be integrated with gamma background)
//double EnergyLossICS(double E, double b)
//{
//	return E*EnergyLossICS(E*b);
//}

/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  /*Inverse Compton Scattering angle averaged differential cross section
    (per e- energy) */
double EMCrosssec::PICSe(double E,double b,double nE)
/* E-cosmic ray electron energy
   b-background photon energy
   nE-outgoing electron energy*/
{
	return PICSph(E, b, E - nE);
}

double EMCrosssec::PICSeR(double Eb, double r)
/*
	Eb=E*b; r = nE/E
	E-cosmic ray electron energy
    b-background photon energy
    nE-outgoing electron energy
*/
{
	return PICSphR(Eb, 1.-r);
}

/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
  /*Pair Production angle averaged differential cross section
   (per e- energy) */
double EMCrosssec::PPPe_old(double E,double b,double nE)
/* E-cosmic ray photon energy
   b-background photon energy
   nE-outgoing electron energy */
{
	if(nE<0) nE=E-fMass;
	double tem1,s1,s2,tem2,tem3;
	tem1=E-nE;            /*e+ energy */
	s1=fSquareMass/nE*E/tem1*E; /*Smin */
	s2=4.0*b*E;
	if(s1>s2) return 0.0;
	tem2=fSquareMass*E*(1/nE+1/tem1);
	tem3=(nE/tem1+tem1/nE)*(s2-s1)+4.0*tem2*log(s2/s1)-4.0*tem2*(1/s1-1/s2)*tem2;
	VERIFY(IsValid(tem3));
	return tem3*C1*8.0/s2/s2/E;
}

/* Pair Production angle averaged differential cross section
//taken from ApJ, 335 p.786 1988, formula (B1)
//Zdziarski article appendix referring to earlier publications of Aharonyan et al.
   E-cosmic ray photon energy
   b-background photon energy
   nE-outgoing electron energy */
double EMCrosssec::PPPeAharonyan(double E,double b,double nE)
{
	double W = E*b/fSquareMass;
	if(W < 1.)
		return 0.;
	double nE2 = E-nE;
	double Wstar = 0.25*E*E/nE/nE2;
	if(Wstar > W || Wstar < 1.)
		return 0.;

	double r = 0.5 * (nE/nE2 + nE2/nE);
	double fracW = Wstar/W;
	double result = r + fracW*(2.*fracW - (2.+r) - 2.*log(fracW));
	ASSERT_VALID_NO(result);
	return (result * 4. * C1 / fSquareMass / E /W);
}

/* Pair Production angle averaged total cross section
//taken from ApJ, 335 p.786 1988, formulas (B11)-(B15)
//Zdziarski article appendix refering to earlier publications of Gould et al.
   E-cosmic ray photon energy
   b-background photon energy
   nE-outgoing electron energy */
double EMCrosssec::RPP_Gould(double E,double b)
{
	double v = E*b/fSquareMass - 1.;
	if(v < 0.)
		return 0.;
	double phi = 0.;
	if(v < 0.01)
	{
		phi = sqrt(v)*v*(1.33333333333333+v*(1.2+v*(-1.72857142857143+v*2.01455026455026)));//formula (B14)
	}
	else if(v > 100)
	{
		double s = 4*E*b/fSquareMass;
		double lns = log(s);
		const double ln4_2 = -0.613705638880109;//log(4)-2
		const double pi_const = 0.289868133696453;//(Pi^2 - 9)/3
		phi = 0.5*s*(lns-2.) + lns*ln4_2 - pi_const + 4./s*(lns + 1.125);//formula (B15)
	}else
	{
		double sqrtV = sqrt(v);
		double sqrtVp1 = sqrt(v+1.);
		double w = (sqrtVp1 + sqrtV)/(sqrtVp1 - sqrtV);// formula (B13)
		double lnw = log(w);
		double lnwp1 = log(w+1.);
		const double pi_const2 = 3.28986813369645;//Pi^2/3
		phi = (1./(1.+v) + 2.*v)*lnw - 2.*sqrtV*(1.+2.*v)/sqrtVp1 - lnw*lnw + 2.*lnwp1*lnwp1 + 4.*li2(1./(1.+w)) - pi_const2;
	}
	ASSERT_VALID_NO(phi);
	double coef = fMass/E/b;
	return (phi * 2. * C1 * coef * coef);
}

/// Integral P(E,b,nE) d(nE)
/// r = nE/E
/// W = E*b/Me2
/// should be multiplied by 4.*C1/Me2/W
double EMCrosssec::intPPPdr(double r, double W)
{
	const double Pi2 = 9.86960440108936;//Pi^2
	double W2 = W*W;
	double lnr = log(r);
	double result = 0;
	if(r < 0.01)
	{
		double ln4W = log(4.*W);
		double W_1 = 1./W;
		double piConst = 1.64493406684823;//Pi^2/6
		result = 0.125/r*W_1*(1.-W_1) + 0.5*lnr;//lowest order term (if W->inf, r->0)
		result += W_1*((0.125*(W_1-1.)-piConst) + 0.5*lnr*(0.5*W_1-1+0.5*lnr+ln4W));//next order term
		result += 0.5*r*(-1 + W_1*(0.75*W_1-3.25+lnr+ln4W));
	}
	else
	{
		double ln1_r = log(1-r);
		double lnLong = 1.38629436111989 + ln1_r + lnr + log(W);// Log(4*(1 - r)*r*W)
		result = ((-1 + W)/(-1 + r) + (-1 + W)/r - 8*r*W2 - 2*(1 - 2*W + 2*W2)*ln1_r + 2*(1 - 2*W + 2*W2)*lnr + 
		 2*W*(-Pi2 + ln1_r*ln1_r - lnr*lnr - 2*ln1_r*lnLong + 2*lnr*lnLong + 2*li2(1. - r) - 2*li2(r)))/(8.*W2);
	}
	return result;
}

/// Integral P(E,b,nE) d(nE) from 0 to r*E
double EMCrosssec::PPPdrIntegral(double E, double b, double r1, double r2)
{
	ASSERT(r1<=r2);
	double W = E*b/fSquareMass;
	if(W<1)
		return 0.;
	double W_1=1/W;
	double rMin = (W_1>1e-2)?0.5*(1.-sqrt(1.-1./W)):0.25*W_1*(1.+0.25*W_1*(1.+0.5*W_1));
	if(r1<=rMin)
		r1 = rMin;
	double rMax = 1.-rMin;
	if(r2>rMax)
		r2 = rMax;
	double result = r1>=r2 ? 0. : (intPPPdr(r2, W)-intPPPdr(r1, W));
	ASSERT_VALID_NO(result);
	return result*4.*C1/fSquareMass*W_1;
}

/// Integral P(E,b,nE) d(nE) from 0 to r*E
double EMCrosssec::lowerPPPdrIntegral(double E, double b, double r)
{
	double W = E*b/fSquareMass;
	if(W<1)
		return 0.;
	double W_1=1/W;
	double rMin = (W_1>1e-2)?0.5*(1.-sqrt(1.-1./W)):0.25*W_1*(1.+0.25*W_1*(1.+0.5*W_1));
	if(r<=rMin)
		return 0.;
	double rMax = 1.-rMin;
	if(r>rMax)
		r = rMax;
	double result = intPPPdr(r, W)-intPPPdr(rMin, W);
	ASSERT_VALID_NO(result);
	return result*4.*C1/fSquareMass*W_1;
}

/// Integral P(E,b,nE) nE d(nE)
/// r = nE/E
/// W = E*b/Me2
/// should be multiplied by 4.*C1/Me2/W*E
double EMCrosssec::intPPPrdr(double r, double W)
{
	const double Pi2 = 9.86960440108936;//Pi^2
	double W2 = W*W;
	double lnr = log(r);
	double result = 0;
	if(r < 0.01)
	{
		double ln4W = log(4.*W);
		result = -(-3. + (3. + 4.*Pi2)*W + 3.*(-1. + W)*lnr)/(24.*W2);
		result += ((1. - 4.*W + 2.*W2 + 2.*W*lnr + 2.*W*ln4W)/(4.*W2)*r);
		result += ((3. - 15.*W - 4.*W2 + 4.*W*lnr + 4.*W*ln4W)/(16.*W2)*r*r);
	}
	else
	{
		double ln1_r = log(1-r);
		double lnLong = 1.38629436111989 + ln1_r + lnr + log(W);// Log(4*(1 - r)*r*W)
		result = ((-1. + W)/(-1. + r) - 4.*r*r*W2 + (-1. + 3.*W - 4.*W2)*ln1_r - 
     (-1. + W)*lnr + 2.*W*(-Pi2 + ln1_r*(ln1_r + 2.*lnr - 2.*lnLong) + 2.*li2(1 - r)
        ))/(8.*W2);
	}
	return result;
}

/// Integral P(E,b,nE) nE d(nE) from 0 to r*E
double EMCrosssec::lowerPPPrdrIntegral(double E, double b, double r)
{
	double W = E*b/fSquareMass;
	if(W<1)
		return 0.;
	double W_1=1/W;
	double rMin = (W_1>1e-2)?0.5*(1.-sqrt(1.-1./W)):0.25*W_1*(1.+0.25*W_1*(1.+0.5*W_1));
	if(r<=rMin)
		return 0.;
	double rMax = 1.-rMin;
	if(r>rMax)
		r = rMax;
	double result = intPPPrdr(r, W)-intPPPrdr(rMin, W);
	ASSERT_VALID_NO(result);
	return result*4.*C1/fSquareMass*W_1*E;
}

double EMCrosssec::PPPe(double E,double b,double nE)
/* E-cosmic ray photon energy
   b-background photon energy
   nE-outgoing electron energy */
{
#ifdef PP_ZDZIARSKI
	return PPPeAharonyan(E, b, nE);
#else

	double r = nE/E;
	double W = E*b;// squared maximal electron energy at CM farme
	if(W<fSquareMass)
		return 0.;//below threshold

	double ratio = fSquareMass/W;

	double beta = 0.;//sqrt(1.-Me2/W);//maximal electron velocity at CM farme
	double beta_1 = 0.;//1.-beta;
	if (ratio<1e-8) {//if Taylor series 1st member is more precise than double
		beta_1 = 0.5*ratio;
		beta = 1.- beta_1;
	}else{
		beta = sqrt(1.-ratio);
		beta_1 = 1.-beta;
	}

	if(r<0.5*beta_1||r>0.5*(1.+beta))
		return 0.;
	double rProd = (1.-r)*r;

	// long expression obtained by mathematica
	double result = 
	((fSquareMass - 4.*rProd*W)*(fSquareMass + (2.*rProd-1.)*W) - 4*fSquareMass*rProd*W*log(fSquareMass/(4*rProd*W)))/(rProd*rProd);

	result *= (0.5*C1/(W*W*W*E));//C1 = Pi/2 * alpha^2 (see above)

	ASSERT_VALID_NO(result);

	//cerr << result << "\t" << PPPeAharonyan(E, b, nE) << endl;

	return result;
#endif
}


/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/

//Series expansion of RICS_Jones for small A
double EMCrosssec::RICS_JonesSeriesZero(double A)
{
	double a0 = 0.333333333333333;//  1/3
	double a1 = -0.222222222222222;// -2/9
	double a2 = 0.2166666666666667;// 13/60
	double a3 = 0.2216666666666667;// 133/600
	double a4 = 0.2269841269841270;// 143/630

	double A2 = A*A;

	return a0 + a1*A + A2*(a2 + a3*A + a4*A2);
}

//Series expansion of RICS_Jones for large A
double EMCrosssec::RICS_JonesSeriesInf(double A)
{
	return (0.5*log(A)-0.25)/A;
}

double EMCrosssec::RICS(double Eb)
/*
	Eb = E*b
	E-cosmic ray electron energy
    b-background photon energy      
    Obtained by integration of PICSph (Johnes formula)
*/
{
	double A = 4.*Eb/fSquareMass;
	double result;

	if(A < 1e-3)
	{
		result =  RICS_JonesSeriesZero(A);
	}
	else if(A > 1e5)
	{
		result =  RICS_JonesSeriesInf(A);
	}
	else
	{
		double Ap = 1.+A;
		result = (-(A*(16. + A*(18. + A))) + 2.*Ap*Ap*(8. + A)*log(Ap) + 
		8.*A*Ap*li2(-A))/(4.*A*A*A*Ap);
	}
	result *= (16.*C1/fSquareMass);
	ASSERT_VALID_NO(result);
	return result;
}

double EMCrosssec::RICS(double E,double b)
/* E-cosmic ray electron energy
   b-background photon energy
*/
{
	return RICS(E * b);
}

    /*Inverse Compton Scattering angle averaged full cross section */
//double EMCrosssec::F_ics(double b, double _m/* 1-b (needed if it is small)*/);

double EMCrosssec::F_I_PICSe(double x,double k)
{
	const double Logx=log(x);
	const double k2=k*k;
	const double x2=x*x;
	const double x3=x2*x;
	const double x1=1./x;
	double result;
	if(x>0.99)//x==1
		result=2.0 - k * ( 21. + 4.* Pi * Pi )/6. - 4. * k * log(k)+8.*(x-1.);
	else 
		result=(2.*k + 2.*k2 - 10.*k*x2 - 2.*k2*x2 + 4.*x3 + k*x3 -
		8.*k*x*log(k/x) - 8.*k*x*log((k*(-1. + x1))/4.)*(-1. + x - Logx)
		+ x*Logx*(8. + 10.*k + 4.*k2 - 8.*k*log(1 - x) + 4.*k*Logx) - 
		8.*k*x*li2(x))/(2*x);
  	return result;
}

double EMCrosssec::I_PICSe(double x_min,double x_max,const double s)
{
//	if(s/Me2<1e-3)
//		return 0.0;//forward scatering must be excluded for stability
	if(x_max>1)
		x_max=1.;
	const double k=fSquareMass/s;
	if(k>1e3)
		return 0.0;//forward scatering must be excluded for stability
	const double x_bound=1./(1+4./k);
	if(x_min<x_bound)
		x_min=x_bound;
	if(x_min>=x_max)
		return 0.;
	double result;
	if(x_min>0.98)
		result=8.*(x_max-x_min);
	else
		result=F_I_PICSe(x_max,k)-F_I_PICSe(x_min,k);
	if(result<0)
	{
		char buf[256];
		sprintf(buf, "I_PICSe<0 : x_min=%.15lg\tx_max=%.15lg\ts=%.15lg\n",x_min,x_max,s);
		WARN(buf);
//		scanf("%*c");
		return 0.;
	}
	return 0.5*C1/s*result;

}

double EMCrosssec::I_PICSph(double x_min,double x_max,const double s)
{
	return I_PICSe(1.-x_max,1.-x_min,s);
}

/*Pair Production angle averaged full cross section */
double EMCrosssec::RPP_old(double E,double b)
/*   E-cosmic ray photon energy
   b-background photon energy
                                  old version?*/
{
	double Li2half=0.582240526465012; /*  li2(1/2)   */
	double tem1,tem2,tem3,tem4;
	tem1=fSquareMass/b/E;
	tem2=sqrt(1-tem1);
	tem4=log(tem1);
	tem3=log(1+tem2)*(-1.5*log(1-tem2)+2.0*ln2+1.0/(1.0-tem2)+1.0/(1.0+tem2)+1.0-tem2*tem2);
	tem3+=-2.0*li2(0.5*(1.0+tem2))+1.0/(1+tem2)+tem2*(0.5*tem2-1.0)+2*Li2half;
	tem3=4.0*fSquareMass*(2.0*tem3-1.0-(2.0*tem4+2.0)/tem1+tem4*tem4+tem1*(1.0-tem4));
	tem3-=4.0*fSquareMass*(2.0*tem2*(1.0/tem1-2.0)+log((2.0*tem2+2.0)/tem1-1.0));
	VERIFY(tem3>=0);
	return 0.5*C1/E/b/E/b*tem3;
}


double EMCrosssec::RPP(double E,double b)
//                        new version (described in diplom)*/
/*   E-cosmic ray photon energy
   b-background photon energy*/
{
#ifdef PP_ZDZIARSKI
	//return RPP_Gould(E, b);
	return lowerPPPdrIntegral(E,b,1.);
#else

	double dif = fSquareMass/E/b;
	if(dif>1.) return 0.;
	double b2=0;//maximum e- velocity in CM frame
	double b_2=0;//1-b2
	double ln_m = 0;//log(1-b2);
	double b4=0;//b2*b2
	if (dif<1e-10){
		b_2 = 0.5*dif;
		b2 = 1.-b_2;
		b4 = 1.-2.*b_2;
		ln_m = log(b_2);
	}else{
		b4=1.-dif;
		b2=sqrt(b4);
		b_2 = 1.-b2;
		ln_m=log(1.-b2);
	}

	double ln_p=log(1+b2);

	double result=4.*(li2(0.5*b_2)-Li2half)-ln_p*ln_p+ln_m*(-4.*ln2+ln_m+2.*ln_p);
	result+=((-2)*b2*(1+b4)+(1+b4*b4)*(ln_p-ln_m))/((b_2)*(1.+b2));
	VERIFY(result>=0);
	return 2.*C1*dif*dif/fSquareMass*result;
#endif
}


