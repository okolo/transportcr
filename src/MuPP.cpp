#include "MuPP.h"
#include "Addfunc.h"
#include <math.h>
#include "Ranges.h"
#include "DecaySpectrum.h"

double CMuPairProduction::PmuPP(double x,double s)
{
	// C1=0.5*Pi/137.035989561/137.035989561=3/16*TompsCS*Me^2=Pi/2*alpha^2
	const double  C1=8.3647043703e-5;

	if((s<Mmu2)||(x<=0.)||(x>=1.))
		return 0.;
	double tmp1=1./x-1.;
	double tmp2=Mmu2*(1./x+1./(1.-x));
	double s1=Mmu2/x/(1.-x); /*Smin */
	double s2=4.0*s;
	if(s1>=s2) return 0.0;
	double result=(tmp1+1./tmp1)*(s2-s1)+4.0*tmp2*log(s2/s1)-4.0*tmp2*tmp2*(1/s1-1/s2);
	VERIFY(IsValid(result));
	return result*C1*8.0/s2/s2;
}

CMuPairProduction::CMuPairProduction(const double bmin, const int mm)
{
	DECLARE_E_SHORTCUT
	//DECLARE_K_SHORTCUT
 const double s = BC().s;
 //extern double PmuPP(double x,double s);
 
 LOG_MESSAGE("Muon pair production coeffitients initialization...");
 MuDecaySpectrum muDecay;
 s_min=(int)(s*log10(Mmu2/bmin/Emin));
 int i,j;
 s_max=mm+nn+1;
 for(j=0;j<2;j++)
 {
 	coef[j]=(double**)Calloc(s_max-s_min, sizeof(void*),"CMuPairProduction");
 	coef[j]-=s_min;
 	for(i=s_min;i<s_max;i++)
 		coef[j][i]=(double*)Calloc(nn, sizeof(double),"CMuPairProduction");
 }	
 
 double S=bmin*Emin*pow(10.,(((double)s_min)+1.)/s);
 const int accurasy=10;
 const double deltaLog=log(10.)/accurasy/s;
 const double mult=pow(10.,1./s/accurasy);
 const double y_max=1/mult;
 for(i=s_min;i<s_max;i++,S*=BC().ss1)
 {
 	double x=Emin/Emax*BC().ss1;
	for(j=0;j<nn-1;j++,x*=BC().ss1)
	{	
		double sum_en=0.;
		double sum_mn=0.;
		for(double y=x*mult;y<=y_max;y*=mult)
		{
			double P = PmuPP(y,S);
			sum_en+=P*muDecay.aen(x/y);
			sum_mn+=P*muDecay.mn(x/y);
		};
		assert(sum_en>=0.);
		assert(sum_mn>=0.);
		coef[en][i][j]=sum_en*deltaLog;
		coef[mn][i][j]=sum_mn*deltaLog;
	}
 }
 LOG_MESSAGE("Done");
}

CMuPairProduction::~CMuPairProduction()
{
 for(int j=0;j<2;j++)
 {
	for(int i=s_min;i<s_max;i++)
 		free(coef[j][i]);
	coef[j]+=s_min;
	free(coef[j]);
 }
}

double CMuPairProduction::GetP(int _nE,int _E,int _b,TProducts product)
{
	const int nn = Ranges().nE();
	int i=_b+_E;
	if(i<s_min)
		return 0.;
	int j=nn-1-(_E-_nE);
	VERIFY((j>=0)&&(j<nn));
	return coef[product][i][j];
}


