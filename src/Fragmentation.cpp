
#include "Fragmentation.h"
#include <math.h>
#include "Addfunc.h"

const double CFragmentation::Nc=3.0;//number of quark colors
const double CFragmentation::Nf=6.0;//number of quark flavors
const double CFragmentation::Leff_MeV=250.0;//
CFragmentation::CFragmentation(double _Ejet/*MeV*/,double a,double b)
{
    double C,Y;
    norm=1.0;
    C=a*a/16.0/b/Nc;
    Y=log(_Ejet/Leff_MeV);
    ksiM=Y*(0.5+sqrt(C/Y)-C/Y);
    sigma2=sqrt(b*Y*Y*Y/36.0/Nc);
    /*double sigma=sqrt(0.5*sigma2);
    norm=1.0/sqrt(2*Pi)/sigma;//Exp(ksiM-0.25*sigma2);*/
	Normalize();
}

double CFragmentation::F(double x)
{
	if(x>1)
		return 0.;
    double Dksi=-1.0*log(x)-ksiM;
    return norm/x*Exp(-Dksi*Dksi/sigma2);
}


CFragmentationQCD::CFragmentationQCD(double _Ejet/*MeV*/):
CFragmentation(_Ejet,103.0/9.0 , 7.0)
{
}

CFragmentationSUSY::CFragmentationSUSY(double _Ejet/*MeV*/):
CFragmentation(_Ejet,11.0 , 3.0)
{
    norm*=0.6;//may be not valid for the case X->q e
}

void CFragmentation::Normalize(double _c)
{
	double x,sum;
	const int ac=1500;
	const double x_min=1e-15;
	double mult=pow(x_min,-0.5/((double)ac));
	double y1=F(x_min)*x_min;
	int i=0;
	for(sum=0.0,x=x_min;i<ac;i++)
	{
		x*=mult;
		double add=4.0*F(x)*x+y1;
		x*=mult;
		y1=F(x)*x;
		add+=y1;
		sum+=add*x;
	};
	sum*=(mult*mult-1.0)/6.0;
	norm*=1.0/sum;
}

