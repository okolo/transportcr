#include "ICS.h"
#include "Ranges.h"
#include "Medium.h"
#include <math.h>
#include "FilePtr.h"
#include "Jacobian.h"
#include "TimeZ.h"
#include "SafeOutput.h"
#include "Function.h"
#include "Units.h"

namespace couplings {

////////////////////  ICSold class members
ICSold::ICSold(void):
	cs(Me)
{
	AddChannel(new Channel_e_e(this, EElectron));
	AddChannel(new Channel_e_e(this, EPositron));
	AddChannel(new Channel_e_gamma(this, EElectron));
	AddChannel(new Channel_e_gamma(this, EPositron));
}

void ICSold::InitCoef(const CBackgroundTable& aBackgr)
{
	CMatrix& be = f_be;
	CMatrix& bph = f_bph;
	CVector& ae = f_a;

	const int mm = aBackgr.nK();
	const int nn = Ranges().nE();
	const double Emin = Ranges().Emin();

	int i,j,k;
	double tem1, tem2, E, E1,sss;
	sss=1.0-BC().ss_1;

	double a1=7.0/90.0;//Bode's integration role
	double a2=16.0/45.0;
	double a3=6.0/45.0;
	double x2=0.5*(1.5/BC().ss1+0.5);
	double x3=0.5*(1.0+BC().ss_1);
	double x4=0.5*(1.5+0.5/BC().ss1);
	double xx2=0.5*(1.5/BC().ss2+0.5);
	double xx3=0.5*(1.0+BC().ss_2);
	double xx4=0.5*(1.5+0.5/BC().ss2);

	for(j=0,E=Emin*BC().ss2;j<nn;j++,E*=BC().ss1)
	{
		for(k=0; k<mm; k++)
		{
			const double Fb=aBackgr.Fmid(k);
			const double b = aBackgr.midBackgroundE()[k];

			tem2=cs.PICSe(E,b,Emin);
			VERIFY(tem2>=0);
			for(E1=Emin*BC().ss1,i=0;i<j;i++,E1*=BC().ss1)
			{
				tem1=cs.PICSe(E,b,E1);
				be[i][j]+=E1*sss*(a1*tem1+a1*tem2+a2*(cs.PICSe(E,b,x2*E1)+cs.PICSe(E,b,x4*E1))+a3*cs.PICSe(E,b,x3*E1))*Fb;
				tem2=tem1;
			};
			tem1=cs.PICSe(E,b,0.99*E);
			be[j][j]+=E*(1.0-BC().ss_2)*(a1*tem1+a1*tem2+a2*(cs.PICSe(E,b,xx2*E)+cs.PICSe(E,b,xx4*E))+a3*cs.PICSe(E,b,xx3*E))*Fb;
		    	// /6.0*(tem2+4.0*cs.PICSe(E,b,E1)+cs.PICSe(E,b,0.99*E))*Fb;
		    	//be[j][j]+=E1*(BC().ss_2-BC().ss_1)*tem2;
			tem2=cs.PICSph(E,b,Emin);
			
			VERIFY(tem2>=0);
			for(E1=Emin*BC().ss1,i=0;i<j;i++,E1*=BC().ss1)
			{
				tem1=cs.PICSph(E,b,E1);
				bph[i][j]+=E1*sss*(a1*tem1+a1*tem2+a2*(cs.PICSph(E,b,x2*E1)+cs.PICSph(E,b,x4*E1))+a3*cs.PICSph(E,b,x3*E1))*Fb;
				tem2=tem1;
			};
			tem1=cs.PICSph(E,b,0.99*E);
			bph[j][j]+=E*(1.0-BC().ss_2)*(a1*tem1+a1*tem2+a2*(cs.PICSph(E,b,xx2*E)+cs.PICSph(E,b,xx4*E))+a3*cs.PICSph(E,b,xx3*E))*Fb;
		}
	}
	for(k=0; k<mm; k++)
	{
		const double Fb=aBackgr.Fmid(k);
		const double b = aBackgr.midBackgroundE()[k];
		for(i=0,E=Emin*BC().ss2;i<nn;i++,E*=BC().ss1)
		{		
			ae[i]+=Fb*cs.RICS(E,b);
			f_intPrdrGamma[i]+=Fb*cs.intPrdrICSph(E, b);
		};
	};
}


void ICSold::AdjustCoef()
{
	CMatrix& be = f_be;
	CMatrix& bph = f_bph;
	CVector& ae = f_a;

	double usageLimit = 0.015;//maximal descrepancy between total sigma and integral of dif sigma for energy conservation based adjustment of coefficients
	const int nn = Ranges().nE();
	int iMin = 1;
	const double binDifRatio = 1.-BC().ss_1;
	int i,j;

	CSafeOutput icsOut;
	if(COEF_TEST_ON)
	{
		icsOut.open(plt_local_dir + "ICSinfoZ" + ToString(redshift.z()));
	}
	icsOut << "E\tEsec\tElead\tratioSec\tratioLead\t<ratioEGamma analytic>\t<ratioEGamma>\t<adjustment status>\n" ;

	for(i=iMin; i<nn; i++)
	{
		double Ei = Ranges().midE()[i];
		ASSERT(ae[i]>=0);
		if(ae[i]==0.)
		{
			icsOut << Ei << "\t no interaction" << endl;
			continue;
		}

		double sumPe = 0;
		double sumEPe = 0;
		double sumPph = 0;
		double sumEPph = 0;
		for(j=0; j<=i; j++)
		{
			double Ej = Ranges().midE()[j];
			ASSERT(be[j][i]>=0.);
			ASSERT(bph[j][i]>=0.);
			sumPe += be[j][i];
			sumPph += bph[j][i];
			sumEPe += Ej*be[j][i];
			sumEPph += Ej*bph[j][i];
		}

		double meanSecElectronEnergy = sumEPe/ae[i];
		double meanSecPhotonEnergy = sumEPph/ae[i];
		double meanAnalyticSecPhotonEnergyRatio = f_intPrdrGamma[i]/ae[i];

		const char* leadingParticle;
		double meanNonLeadingEnergy, meanLeadingEnergy, sumPsec, sumEPsec, sumPlead, sumEPlead;
		CMatrix* pbLeading = 0;
		CMatrix* pbSecondary = 0;
		double mostAcurateNonAnaliticalRsecGamma;
		if(meanAnalyticSecPhotonEnergyRatio<0.5)
		//if(meanSecElectronEnergy>meanSecPhotonEnergy)
		{
			leadingParticle = "electron";
			meanNonLeadingEnergy = meanSecPhotonEnergy;
			meanLeadingEnergy = meanSecElectronEnergy;
			pbLeading = &be;
			pbSecondary = &bph;
			sumPsec = sumPph;
			sumEPsec = sumEPph;
			sumPlead = sumPe;
			sumEPlead = sumEPe;
			mostAcurateNonAnaliticalRsecGamma = meanSecPhotonEnergy/Ei;
		}
		else
		{
			leadingParticle = "photon";
			meanNonLeadingEnergy = meanSecElectronEnergy;
			meanLeadingEnergy = meanSecPhotonEnergy;
			pbLeading = &bph;
			pbSecondary = &be;
			sumPsec = sumPe;
			sumEPsec = sumEPe;
			sumPlead = sumPph;
			sumEPlead = sumEPph;
			mostAcurateNonAnaliticalRsecGamma = 1.-meanSecElectronEnergy/Ei;
		}
		CMatrix& bLeading = *pbLeading;
		CMatrix& bSecondary = *pbSecondary;
		double ratioSec = sumPsec/ae[i];//!if ratio > 1 it may mean error
		double ratioLead = sumPlead/ae[i];
		

		icsOut << Ei << "\t" << meanNonLeadingEnergy << "\t" << meanLeadingEnergy <<"\t" << ratioSec << "\t" << ratioLead << "\t"
			<< meanAnalyticSecPhotonEnergyRatio << "\t" << mostAcurateNonAnaliticalRsecGamma <<"\t"<< leadingParticle << "Leading";
		
		if(meanAnalyticSecPhotonEnergyRatio<binDifRatio)//most of leading electrons are concentrated inside the highest bin (low energy limit)
		{
			sumPe -= (be[i][i] + be[i-1][i]);
			sumEPe -= Ei*(be[i][i] + BC().ss_1*be[i-1][i]);
			double bi_1_i = (Ei*(ae[i]*meanAnalyticSecPhotonEnergyRatio-sumPe) + sumEPe)/(Ei*binDifRatio);
			sumPe += bi_1_i;
			double bi_i = ae[i]-sumPe;
			if(bi_1_i>=0&&bi_i>=0)
			{
				be[i-1][i] = bi_1_i;
				be[i][i] = bi_i;
				icsOut << "Ok:SmallEsecMode";
			}else
			{
				ASSERT(false);
				icsOut << "Failed:SmallEsecMode!!!";
			}
		}
		else if(1.-ratioSec<usageLimit)//enough information for energy conservation based adjustment
		{	
			double multiplier = 1./ratioSec;
			for(j=0; j<=i; j++)
			{
				bSecondary[j][i] *= multiplier;
			}
			sumEPsec *= multiplier;
			int k;
			for(k=i;bLeading[i][i]==0.;k--);
			double Ek = Ranges().midE()[k];
			sumPlead -= (bLeading[k][i] + bLeading[k-1][i]);
			sumEPlead -= Ek*(bLeading[k][i] + BC().ss_1*bLeading[k-1][i]);
			double bk_1_i = (ae[i]*(Ek-Ei) + sumEPsec+sumEPlead-Ek*sumPlead)/(Ek-Ranges().midE()[k-1]);
			ASSERT(bk_1_i>=0);
			sumPlead += bk_1_i;
			double bk_i = ae[i] - sumPlead;
			
			if(bk_1_i>=0&&bk_i>=0)
			{
				bLeading[k-1][i] = bk_1_i;
				bLeading[k][i] = bk_i;
				sumEPlead += Ek*(bLeading[k][i] + BC().ss_1*bLeading[k-1][i]);
				icsOut << "Ok: rSec=" << sumEPsec/ae[i]/Ei << " rLead="<<sumEPlead/ae[i]/Ei<<" i-k=" << (i-k);
			}else
			{
				ASSERT(false);
				icsOut << "Failed i-k=" << (i-k);
			}
		}
		else
		{
			icsOut << "Skipped";
		}
		icsOut << endl;
	};
	icsOut.close();
}

void ICSold::Channel_e_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CMatrix& be = fCoupling.f_be;
	const CVector& ae = fCoupling.f_a;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		if(ae[iPrim]>0)
		{// use dif cross section
			aCoef.Add(iPrim, iPrim, -ae[iPrim]);
			for(int iSec = 0; iSec<=iPrim; iSec++)
			{
				aCoef.Add(iSec, iPrim, be[iSec][iPrim]);
			}
		}
	}
}

void ICSold::Channel_e_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const CMatrix& bph = fCoupling.f_bph;
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(int iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, bph[iSec][iPrim]);
	}
}



void ICSold::SetBackgrounds(const Medium& aPropagCoef)
{
	const bool LEE_ADJUSTMENT_ICS = false;//coef adjustments (experimental)
	const int nn = Ranges().nE();
	Coupling::SetBackgrounds(aPropagCoef);
	if(!f_a.length())
	{
		f_intPrdrGamma.create(nn);
		f_a.create(nn);
		f_be.create(nn);
		f_bph.create(nn);
	}
	else
	{
		f_intPrdrGamma.reset();
		f_a.reset();
		f_be.reset();
		f_bph.reset();
	}
	InitCoef(*(aPropagCoef.background()));
	if(LEE_ADJUSTMENT_ICS)
		AdjustCoef();
	FinalCoefAdjustment();
}

void ICSold::DebugOutput(const char* aFolder) const
{
	ASSERT(f_a.length()>0);//coefficients must be initialized
	const int nn = Ranges().nE();
	double ycoeff=units.Lunit/units.Mpc_cm;
	double E = Ranges().Emin()*units.Eunit/units.outEunitMeV;
	CFilePtr eInt(Fopen("icsIntLengthC", aFolder, "wt"));
	for(int i=0; i<nn; i++,E*=BC().ss1)
	{
		if(f_a[i]>0)
				fprintf(eInt,"%lg\t%lg\n", E, ycoeff/f_a[i]);
	}
}

void ICSold::FinalCoefAdjustment()
{
	const int nn = Ranges().nE();
	for(int j=0;j<nn;j++)
	{
		double sum1=0;
		for(int i=0;i<=j;i++)
		{
			sum1+=f_be[i][j];
		};
		if(sum1>f_a[j])
		{
			if((sum1-f_a[j])/sum1>1e-10)
			{
				ASSERT(false);
				LOG_ERROR_ONCE("ae[j]<sum(be[i][j])");
			}
			f_a[j]=sum1;
		}
	}
}

////////////////////  ICScel class members

ICScel::ICScel(void):
	cs(Me)
{
	READ_DOUBLE_SETTING(ICS_CEL_MAX_BINS_DIFF);

	AddChannel(new Channel_e_e(this, EElectron));
	AddChannel(new Channel_e_e(this, EPositron));
	AddChannel(new Channel_e_gamma(this, EElectron));
	AddChannel(new Channel_e_gamma(this, EPositron));
}


double ICScel::PICSe(double E,double b,double nE)
{
	EMCrosssec cs(ParticleData::getParticleMass(EElectron));
	return cs.PICSe(E,b,nE);
}

double ICScel::PICSph(double E,double b,double nE)
{
	EMCrosssec cs(ParticleData::getParticleMass(EElectron));
	return cs.PICSph(E,b,nE);
}

void ICScel::InitCoef(const CBackgroundTable& aBackgr)
{
	const int nn = Ranges().nE();
	const int mm = aBackgr.nK();
	const CBinning midE = Ranges().midE();
	const CBinning& midK = aBackgr.midBackgroundE();
	const double Emin = Ranges().Emin();
	{
		int jMax = 2*nn + 1;
		double E=Emin;
		for(int j=0; j<jMax; j++, E*=BC().ss2)
		{
			double sum = 0;
			for(int k=0; k<mm; k++)
				sum += aBackgr.Fmid(k)*cs.EnergyLossICS(E*midK[k]);
			fEnergyLossRate[j] = sum;
		}
	}

	CMatrix& be = f_be;
	CMatrix& bph = f_bph;
	CVector& ae = f_a;

	PFunction func;
	double maxAbsError=0;
	double maxRelError=1e-2;
	size_t maxSubIntervals=1000;
	for(int k=0; k<mm; k++)
	{
		double Fb=aBackgr.Fmid(k);
		double b = midK[k];
		func.b = b;
		double E=Emin*BC().ss2;
		for(int j=0;j<nn;j++,E*=BC().ss1)
		{
			func.E = E;
			double aK=Fb*cs.RICS(E,b);
			ae[j]+=aK;
			f_intPrdrGamma[j]+=Fb*cs.intPrdrICSph(E, b);
#ifdef ICS_DEBUG
			double eLossK = Fb*cs.EnergyLossICS(E*b);
			double deltaErel = eLossK/aK;//typical fraction of energy lost in single interaction
#endif
			{//be
				func.P = PICSe;
#ifdef ICS_DEBUG
				double deltaEbins = -log(1.-deltaErel)/log(BC().ss1);//mean difference in bins between primary and secondary electron
				double Psum = 0;
				double PEsum = 0;
				int iLast = -1;
				int iFirst = -1;
#endif
				double E1=Emin*BC().ss1;
				for(int i=0;i<=j;i++,E1*=BC().ss1)
				{
					double bK = Fb*funcUtils.gsl_integration_qag(&func, E1*BC().ss_1, E1<E ? E1 : E, maxAbsError, maxRelError, maxSubIntervals); //E1*sss*(a1*tem1+a1*tem2+a2*(cs.PICSe(E,b,x2*E1)+cs.PICSe(E,b,x4*E1))+a3*cs.PICSe(E,b,x3*E1))*Fb;
					ASSERT_VALID_NO(bK);
					be[i][j] += bK;
#ifdef ICS_DEBUG
					if(bK>0)
					{
						iLast = i;
						if(iFirst<0)
							iFirst = i;
						Psum += bK;
						PEsum += bK*midE[i];
					}
#endif
				};
#ifdef ICS_DEBUG
				if(j>2*deltaEbins && iLast>0 && iFirst<iLast)
				{
					double PEsumTh = midE[j]*(aK-eLossK);
					double ratioA = aK/Psum;
					double ratioPE = PEsumTh/PEsum;
					ASSERT(fabs(ratioA-1.)<5e-2);
					ASSERT(fabs(ratioPE-1.)<5e-2);
				}
#endif
			}
			{//bph
				func.P = PICSph;
#ifdef ICS_DEBUG
				double deltaEbins = -log(deltaErel)/log(BC().ss1);//mean difference in bins between primary electron and secondary photon
				double Psum = 0;
				double PEsum = 0;
				int iLast = -1;
				int iFirst = -1;
#endif
				double E1=Emin*BC().ss1;
				for(int i=0;i<=j;i++,E1*=BC().ss1)
				{
					double bK = Fb*funcUtils.gsl_integration_qag(&func, E1*BC().ss_1, E1<E ? E1 : E, maxAbsError, maxRelError, maxSubIntervals); //E1*sss*(a1*tem1+a1*tem2+a2*(cs.PICSe(E,b,x2*E1)+cs.PICSe(E,b,x4*E1))+a3*cs.PICSe(E,b,x3*E1))*Fb;
					ASSERT_VALID_NO(bK);
					bph[i][j] += bK;
#ifdef ICS_DEBUG
					if(bK>0)
					{
						iLast = i;
						if(iFirst<0)
							iFirst = i;
						Psum += bK;
						PEsum += bK*midE[i];
					}
#endif
				};
#ifdef ICS_DEBUG
				if(j>4*deltaEbins && iLast>0 && iFirst<iLast)
				{
					double PEsumTh = midE[j]*eLossK;
					double ratioA = aK/Psum;
					double ratioPE = PEsumTh/PEsum;
					ASSERT(fabs(ratioA-1.)<5e-2);
					ASSERT(fabs(ratioPE-1.)<5e-2);
				}
#endif
			}
		}
	}

#ifdef ICS_DEBUG
	{
		double ycoeff=units.Lunit/units.Mpc;
		double E;
		int i;
		
		CFilePtr eInt(Fopen("icsIntLengthOrig", plt_local_c, "wt"));
		for(i=0,E=Emin*units.Eunit/units.outEunitMeV;i<nn;i++,E*=BC().ss1)
		{
			if(ae[i]>0)
					fprintf(eInt,"%lg\t%lg\n", E, ycoeff/ae[i]);
		}
	}
#endif

}

double ICScel::ICS_CEL_MAX_BINS_DIFF=0.;//maximal difference in bins between final and initial energy of electron in ICS for which CEL approximation is used
//if ICS_CEL_MAX_BINS_DIFF=0 optimal value of ICS_CEL_MAX_BINS_DIFF is found based on effective interaction rate minimization
//if ICS_CEL_MAX_BINS_DIFF<0 CEL approximation is not used but coefficient are calculated by ICScel class

void ICScel::AdjustCoef()
{
	CVector& ae = f_a;
	double usageLimit = 0.015;//maximal discrepancy between total sigma and integral of dif sigma for energy conservation based adjustment of coefficients

	const int nn = Ranges().nE();
	const CBinning& midE = Ranges().midE();
	int iMin = 1;
	const double binDifRatio = 1.-BC().ss_1;
	//int i,j;

	if(ICS_CEL_MAX_BINS_DIFF>=0.)
	{
		{//f_ph adjustment
			CMatrix& b = f_bph;
			for(int i=iMin; i<nn; i++)
			{
				if(ae[i]>0.)
				{
					double eLoss = fEnergyLossRate[2*i+1];
					double deltaErel = eLoss/ae[i];
					double deltaEbins = -log(deltaErel)/log(BC().ss1);//difference in bins between energy of primary electron and secondary photon
					if(i<4*deltaEbins)
						continue;//not enough coefs for adjustment
					int jMeanSec = i-(int)deltaEbins;
					double Emean = midE[jMeanSec];
					CMatrix system(2);
					double PEsumTh = midE[i]*eLoss;

					int j=0;
					for(; j<=jMeanSec; j++)
					{
						double P = b[j][i];
						system[0][0] += P*Emean;
						system[1][0] += P*midE[j];
					}
					for(; j<=i; j++)
					{
						double P = b[j][i];
						system[0][1] += P*Emean;
						system[1][1] += P*midE[j];
					}
#ifdef ICS_DEBUG
					double ratioAorig = ae[i]*Emean/(system[0][0]+system[0][1]);//needed just for debugging
					double ratioPEorig = PEsumTh/(system[1][0]+system[1][1]);//needed just for debugging
					ASSERT(fabs(ratioAorig-1.)<5e-2);
					ASSERT(fabs(ratioPEorig-1.)<5e-2);
#endif
					CVector x(2);
					x[0] = ae[i]*Emean;
					x[1] = PEsumTh;
					system.solveLinearSystem(x);
					for(int l=0; l<=1; l++)
						ASSERT(fabs(x[l]-1.)<5e-2);
					for(j=0;j<=jMeanSec;j++)
						b[j][i] *= x[0];
					for(;j<=i;j++)
						b[j][i] *= x[1];
#ifdef ICS_DEBUG
					double Psum = system[0][0]*x[0] + system[0][1]*x[1];
					double PEsum = system[1][0]*x[0] + system[1][1]*x[1];
					double ratioA = ae[i]*Emean/Psum;
					double ratioPE = PEsumTh/PEsum;
					ASSERT(fabs(ratioA-1.)<1e-5);
					ASSERT(fabs(ratioPE-1.)<1e-5);
#endif
				}
			}
		}

		{//f_be adjustment
			CMatrix& b = f_be;
			for(int i=iMin; i<nn; i++)
			{
				if(ae[i]>0.)
				{
					double eLoss = fEnergyLossRate[2*i+1];
					double deltaErel = eLoss/ae[i];
					double deltaEbins = -log(1.-deltaErel)/log(BC().ss1);//difference in bins between energy of primary and secondary electron
					if(i<4*deltaEbins)
						continue;//not enough coefs for adjustment

					bool UseCEL = ICS_CEL_MAX_BINS_DIFF>=0;
					if(UseCEL)
					{
						if(ICS_CEL_MAX_BINS_DIFF==0.)//using CEL approximation only if it minimizes effective absorption coef
						{
							double coef = Jacobian::GetCELintRate(fEnergyLossRate, i);
							UseCEL = (coef < ae[i]);
						}
						else
						{//use CEL if difference in bins is smaller then ICS_CEL_BINS_MAX_DIFF
							UseCEL = (deltaEbins < ICS_CEL_MAX_BINS_DIFF);
						}
					}
					if(UseCEL)
						continue;//no need to adjust coefficients, since they are not being used

					int jMeanSec = i-(int)deltaEbins;
					double Emean = midE[jMeanSec];
					CMatrix system(2);
					double PEsumTh = midE[i]*(ae[i]-eLoss);

					int j=0;
					for(; j<=jMeanSec; j++)
					{
						double P = b[j][i];
						system[0][0] += P*Emean;
						system[1][0] += P*midE[j];
					}
					for(; j<=i; j++)
					{
						double P = b[j][i];
						system[0][1] += P*Emean;
						system[1][1] += P*midE[j];
					}


#ifdef ICS_DEBUG
					double ratioAorig = ae[i]*Emean/(system[0][0]+system[0][1]);//needed just for debugging
					double ratioPEorig = PEsumTh/(system[1][0]+system[1][1]);//needed just for debugging
					ASSERT(fabs(ratioAorig-1.)<0.1);//PROBLEM: should be smaller
					ASSERT(fabs(ratioPEorig-1.)<0.2);//PROBLEM: should be smaller
#endif
					CVector x(2);
					x[0] = ae[i]*Emean;
					x[1] = PEsumTh;
					system.solveLinearSystem(x);

#ifdef ICS_DEBUG
					for(int l=0; l<=1; l++)
						ASSERT(fabs(x[l]-1.)<0.3);//PROBLEM: should be smaller (also see values of fabs(ratioAorig-1.) fabs(ratioPEorig-1.) above)
#pragma message ("Lee Adjustment problem if Sarkar background is used")
#endif
					for(j=0;j<=jMeanSec;j++)
						b[j][i] *= x[0];
					for(;j<=i;j++)
						b[j][i] *= x[1];
#ifdef ICS_DEBUG
					double Psum = system[0][0]*x[0] + system[0][1]*x[1];
					double PEsum = system[1][0]*x[0] + system[1][1]*x[1];
					double ratioA = ae[i]*Emean/Psum;
					double ratioPE = PEsumTh/PEsum;
					ASSERT(fabs(ratioA-1.)<1e-5);
					ASSERT(fabs(ratioPE-1.)<1e-5);
#endif
				}
			}
		}
		return;
	}/**/

	CMatrix& be = f_be;
	CMatrix& bph = f_bph;

	CSafeOutput icsOut;
	if(COEF_TEST_ON)
	{
		icsOut.open(plt_local_dir + "ICSinfoZ" + ToString(redshift.z()));
	}
	icsOut << "E\tEsec\tElead\tratioSec\tratioLead\t<ratioEGamma analytic>\t<ratioEGamma>\t<adjustment status>\n" ;

	for(int i=iMin; i<nn; i++)
	{
		double Ei = Ranges().midE()[i];
		ASSERT(ae[i]>=0);
		if(ae[i]==0.)
		{
			icsOut << Ei << "\t no interaction" << endl;
			continue;
		}

		double sumPe = 0;
		double sumEPe = 0;
		double sumPph = 0;
		double sumEPph = 0;
		for(int j=0; j<=i; j++)
		{
			double Ej = Ranges().midE()[j];
			ASSERT(be[j][i]>=0.);
			ASSERT(bph[j][i]>=0.);
			sumPe += be[j][i];
			sumPph += bph[j][i];
			sumEPe += Ej*be[j][i];
			sumEPph += Ej*bph[j][i];
		}

		double meanSecElectronEnergy = sumEPe/ae[i];
		double meanSecPhotonEnergy = sumEPph/ae[i];
		double meanAnalyticSecPhotonEnergyRatio = f_intPrdrGamma[i]/ae[i];

		const char* leadingParticle;
		double meanNonLeadingEnergy, meanLeadingEnergy, sumPsec, sumEPsec, sumPlead, sumEPlead;
		CMatrix& bLeading = be;
		CMatrix& bSecondary = bph;
		double mostAcurateNonAnaliticalRsecGamma;
		if(meanAnalyticSecPhotonEnergyRatio<0.5)
		//if(meanSecElectronEnergy>meanSecPhotonEnergy)
		{
			leadingParticle = "electron";
			meanNonLeadingEnergy = meanSecPhotonEnergy;
			meanLeadingEnergy = meanSecElectronEnergy;
			sumPsec = sumPph;
			sumEPsec = sumEPph;
			sumPlead = sumPe;
			sumEPlead = sumEPe;
			mostAcurateNonAnaliticalRsecGamma = meanSecPhotonEnergy/Ei;
		}
		else
		{
			leadingParticle = "photon";
			meanNonLeadingEnergy = meanSecElectronEnergy;
			meanLeadingEnergy = meanSecPhotonEnergy;
			bLeading = bph;
			bSecondary = be;
			sumPsec = sumPe;
			sumEPsec = sumEPe;
			sumPlead = sumPph;
			sumEPlead = sumEPph;
			mostAcurateNonAnaliticalRsecGamma = 1.-meanSecElectronEnergy/Ei;
		}
		double ratioSec = sumPsec/ae[i];//!if ratio > 1 it may mean error
		double ratioLead = sumPlead/ae[i];
		

		icsOut << Ei << "\t" << meanNonLeadingEnergy << "\t" << meanLeadingEnergy <<"\t" << ratioSec << "\t" << ratioLead << "\t"
			<< meanAnalyticSecPhotonEnergyRatio << "\t" << mostAcurateNonAnaliticalRsecGamma <<"\t"<< leadingParticle << "Leading";
		
		if(meanAnalyticSecPhotonEnergyRatio<binDifRatio)//most of leading electrons are concentrated inside the highest bin (low energy limit)
		{
			sumPe -= (be[i][i] + be[i-1][i]);
			sumEPe -= Ei*(be[i][i] + BC().ss_1*be[i-1][i]);
			double bi_1_i = (Ei*(ae[i]*meanAnalyticSecPhotonEnergyRatio-sumPe) + sumEPe)/(Ei*binDifRatio);
			sumPe += bi_1_i;
			double bi_i = ae[i]-sumPe;
			if(bi_1_i>=0&&bi_i>=0)
			{
				be[i-1][i] = bi_1_i;
				be[i][i] = bi_i;
				icsOut << "Ok:SmallEsecMode";
			}else
			{
				ASSERT(false);
				icsOut << "Failed:SmallEsecMode!!!";
			}
		}
		else if(1.-ratioSec<usageLimit)//enough information for energy conservation based adjustment
		{	
			double multiplier = 1./ratioSec;
			for(int j=0; j<=i; j++)
			{
				bSecondary[j][i] *= multiplier;
			}
			sumEPsec *= multiplier;
			int k;
			for(k=i;bLeading[i][i]==0.;k--);
			double Ek = Ranges().midE()[k];
			sumPlead -= (bLeading[k][i] + bLeading[k-1][i]);
			sumEPlead -= Ek*(bLeading[k][i] + BC().ss_1*bLeading[k-1][i]);
			double bk_1_i = (ae[i]*(Ek-Ei) + sumEPsec+sumEPlead-Ek*sumPlead)/(Ek-Ranges().midE()[k-1]);
			ASSERT(bk_1_i>=0);
			sumPlead += bk_1_i;
			double bk_i = ae[i] - sumPlead;
			
			if(bk_1_i>=0&&bk_i>=0)
			{
				bLeading[k-1][i] = bk_1_i;
				bLeading[k][i] = bk_i;
				sumEPlead += Ek*(bLeading[k][i] + BC().ss_1*bLeading[k-1][i]);
				icsOut << "Ok: rSec=" << sumEPsec/ae[i]/Ei << " rLead="<<sumEPlead/ae[i]/Ei<<" i-k=" << (i-k);
			}else
			{
				ASSERT(false);
				icsOut << "Failed i-k=" << (i-k);
			}
		}
		else
		{
			icsOut << "Skipped";
		}
		icsOut << endl;
	};
	icsOut.close();
}


void ICScel::Channel_e_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CVector& eLoss = fCoupling.fEnergyLossRate;
	const CMatrix& be = fCoupling.f_be;
	const CVector& ae = fCoupling.f_a;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		if(ae[iPrim]<=0)
			continue;
		//double E = Ranges().midE()[iPrim];
		double coef = Jacobian::GetCELintRate(eLoss, iPrim);
		bool UseCEL = ICS_CEL_MAX_BINS_DIFF>=0;
		if(UseCEL)
		{
			if(ICS_CEL_MAX_BINS_DIFF==0.)//using CEL approximation only if it minimizes effective absorption coef
				UseCEL = (coef < ae[iPrim]);
			else
			{//use CEL if difference in bins is smaller then ICS_CEL_BINS_MAX_DIFF
				double deltaErel = eLoss[2*iPrim+1]/ae[iPrim];//typical energy lost in single interaction
				double deltaEbins = log(1/(1-deltaErel))/log(BC().ss1);//difference in bins
				UseCEL = deltaEbins < ICS_CEL_MAX_BINS_DIFF;
			}
		}
		if(UseCEL)
		{//use CEL approximation
			Jacobian::AddCELbin(aCoef, eLoss, iPrim);
		}
		else
		{// use dif cross section
			aCoef.Add(iPrim, iPrim, -ae[iPrim]);
			for(int iSec = 0; iSec<=iPrim; iSec++)
			{
				aCoef.Add(iSec, iPrim, be[iSec][iPrim]);
			}
		}
	}
}

void ICScel::Channel_e_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const CMatrix& bph = fCoupling.f_bph;
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(int iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, bph[iSec][iPrim]);
	}
}



void ICScel::SetBackgrounds(const Medium& aPropagCoef)
{
	const bool LEE_ADJUSTMENT_ICS = false;//coef adjustments (experimental)
	Coupling::SetBackgrounds(aPropagCoef);

	const int nn = Ranges().nE();
	int jMax = 2*nn + 1;

	if(!fEnergyLossRate.length())
	{
		fEnergyLossRate.create(jMax);
		f_intPrdrGamma.create(nn);
		f_a.create(nn);
		f_be.create(nn);
		f_bph.create(nn);
	}
	else
	{
		fEnergyLossRate.reset();
		f_intPrdrGamma.reset();
		f_a.reset();
		f_be.reset();
		f_bph.reset();
	}

	InitCoef(*(aPropagCoef.background()));

	if(LEE_ADJUSTMENT_ICS)
		AdjustCoef();
	//FinalCoefAdjustment();

	double E=Ranges().Emin();
	if(COEF_TEST_ON && redshift.z()==0.)
	{
		CFilePtr file(FopenL(plt_local_dir + "ICS_CEL","wt"));
		for(int j=0; j<jMax; j++, E*=BC().ss2)
		{//rate in Mpc^-1
			fprintf(file, "%lg\t%lg\n", (E*units.Eunit*1e6), (units.Mpc_cm /units.Lunit*fEnergyLossRate[j]));
		}
	}
}

void ICScel::DebugOutput(const char* aFolder) const
{
	ASSERT(f_a.length()>0);//coefficients must be initialized
	const int nn = Ranges().nE();
	double ycoeff=units.Lunit/units.Mpc_cm;
	double E = Ranges().Emin()*units.Eunit/units.outEunitMeV;
	CFilePtr eInt(Fopen("icsIntLengthC", aFolder, "wt"));
	for(int i=0; i<nn; i++,E*=BC().ss1)
	{
		if(f_a[i]>0)
				fprintf(eInt,"%lg\t%lg\n", E, ycoeff/f_a[i]);
	}
	if(ICS_CEL_MAX_BINS_DIFF>=0.)
	{
		CFilePtr eL(Fopen("icsElossLengthC", aFolder, "wt"));
		int jMax = 2*nn+1;
		E=Ranges().Emin()*units.Eunit/units.outEunitMeV;
		for(int j=0; j<jMax; j++, E*=BC().ss2)
		{
			if(fEnergyLossRate[j]>0)
				fprintf(eL,"%lg\t%lg\n", E, ycoeff/fEnergyLossRate[j]);
		}
	}
}

////////////////////  ICSsplitted class members

double ICSsplitted::ICS_CEL_SPLITTED_MAX_BINS_DIFF=1.;//maximal difference in bins between final and initial energy of electron in ICS for which CEL approximation is used

ICSsplitted::ICSsplitted():
fRcel(0),
fRdiscr(0),
fRho(0),
cs(Me)
{
	//READ_DOUBLE_SETTING(ICS_CEL_SPLITTED_MAX_BINS_DIFF);
	AddChannel(new Channel_e_e(this, EElectron));
	AddChannel(new Channel_e_e(this, EPositron));
	AddChannel(new Channel_e_gamma(this, EElectron));
	AddChannel(new Channel_e_gamma(this, EPositron));
}

void ICSsplitted::Channel_e_e::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	const CVector& eLoss = fCoupling.fEnergyLossRate;
	const CMatrix& be = fCoupling.fb_discr;
	const CVector& ae = fCoupling.fa_discr;

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		//CEL part
		Jacobian::AddCELbin(aCoef, eLoss, iPrim);

		//Discrete part
		if(ae[iPrim]>0)
		{
			aCoef.Add(iPrim, iPrim, -ae[iPrim]);
			for(int iSec = 0; iSec<=iPrim; iSec++)
			{
				aCoef.Add(iSec, iPrim, be[iSec][iPrim]);
			}
		}
	}
}

void ICSsplitted::Channel_e_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
		for(int iSec = 0; iSec<=iPrim; iSec++)
			aCoef.Add(iSec, iPrim, fCoupling.f_bph[iSec][iPrim]);
	}
}


ICSsplitted::~ICSsplitted()
{
	delete fRcel;
	delete fRdiscr;
	delete fRho;
}

double ICSsplitted::PICSphR(double Eb, double r)
{
	EMCrosssec cs(ParticleData::getParticleMass(EElectron));
	return cs.PICSphR(Eb, r);
}

double ICSsplitted::PICSeR(double Eb, double r)
{
	EMCrosssec cs(ParticleData::getParticleMass(EElectron));
	return cs.PICSeR(Eb, r);
}

void ICSsplitted::Init(double aKmin, double aKmax)
{
	const int accuracy = 10000;
	double wMin = Ranges().Emin()*aKmin;
	double w = wMin;
	double wMax=Ranges().Emax()*aKmax;
	double rMinCEL = pow(10.,-ICS_CEL_SPLITTED_MAX_BINS_DIFF/BC().s);
	CVarVector Rcel;
	CVarVector Rdiscr;
	CVarVector W;
	CVarVector Rho;//mean Eph/E for CEL part
	for(int outputCounter=0; w<=wMax;w*=BC().ss2, outputCounter++)
	{
		W.add(w);
		RKern P(PICSphR, w);
		double Rtot = cs.RICS(w);
		double rMin = 1.-cs.rMaxICSph(w);
		double celR, discrR, rho;
		if(rMin < rMinCEL)
		{
			celR=FuncUtils::integrate(0, 1.-rMinCEL, &P, &gUnitFunc, accuracy);
			ASSERT_VALID_NO(celR);
			discrR = Rtot-celR;
			ASSERT_VALID_NO(discrR);
			rho = celR>0 ? FuncUtils::integrate(0., 1.-rMinCEL, &P, &gLinearFunc, accuracy)/celR : 0.;
			ASSERT_VALID_NO(rho);
		}
		else
		{
			celR = Rtot;
			discrR = 0.;
			rho = cs.EnergyLossICS(w)/Rtot;
			ASSERT_VALID_NO(rho);
		}
		Rcel.add(celR);
		Rdiscr.add(discrR);
		Rho.add(rho);
//		if(outputCounter%10 == 0)
//			cout << 4.*w/Me2 << '\t' << celR/Rtot << '\t' << Rtot << '\t' << rho << endl;
	}
	fRcelV.create(Rcel);
	fRdiscrV.create(Rdiscr);
	fW.create(W);
	fRhoV.create(Rho);
	fRcel = new CLinearFunc(fW, fRcelV);
	fRdiscr = new CLinearFunc(fW, fRdiscrV);
	fRho = new CLinearFunc(fW, fRhoV);
	const int nn = Ranges().nE();
	fa_discr.create(nn);
	fb_discr.create(nn, nn);
	f_bph.create(nn);
}

void ICSsplitted::InitCoef(const CBackgroundTable& aBackgr)
{
	const double Emin = Ranges().Emin();
	const int nn = Ranges().nE();

	const double bmin = aBackgr.Kmin();
	const int mm = aBackgr.nK();
	int i,j,k;

	CMatrix& bph = f_bph;
	bph.reset();

	double tem1, tem2, E, E1, b, sss, Fb;
	sss=1.0-BC().ss_1;

	double a1=7.0/90.0;//Bode's integration role
	double a2=16.0/45.0;
	double a3=6.0/45.0;
	double x2=0.5*(1.5/BC().ss1+0.5);
	double x3=0.5*(1.0+BC().ss_1);
	double x4=0.5*(1.5+0.5/BC().ss1);
	double xx2=0.5*(1.5/BC().ss2+0.5);
	double xx3=0.5*(1.0+BC().ss_2);
	double xx4=0.5*(1.5+0.5/BC().ss2);

	for(j=0,E=Emin*BC().ss2;j<nn;j++,E*=BC().ss1)
	{
		for(k=0,b=bmin*BC().ss2;k<mm;k++,b*=BC().ss1)
		{
			Fb=aBackgr.Fmid(k);
			tem2=cs.PICSph(E,b,Emin);
			VERIFY(tem2>=0);
			for(E1=Emin*BC().ss1,i=0;i<j;i++,E1*=BC().ss1)
			{
				tem1=cs.PICSph(E,b,E1);
				bph[i][j]+=E1*sss*(a1*tem1+a1*tem2+a2*(cs.PICSph(E,b,x2*E1)+cs.PICSph(E,b,x4*E1))+a3*cs.PICSph(E,b,x3*E1))*Fb;
				tem2=tem1;
			};
			tem1=0.; // = cs.PICSph(E,b,E) // originaly it was replaced by cs.PICSph(E,b,0.99*E), which is wrong
			bph[j][j]+=E*(1.0-BC().ss_2)*(a1*tem1+a1*tem2+a2*(cs.PICSph(E,b,xx2*E)+cs.PICSph(E,b,xx4*E))+a3*cs.PICSph(E,b,xx3*E))*Fb;
		}
	}
}

void ICSsplitted::SetBackgrounds(const Medium& aPropagCoef)
{
	Coupling::SetBackgrounds(aPropagCoef);
	ASSERT(ICS_CEL_SPLITTED_MAX_BINS_DIFF>0.);
	const CBackgroundTable& backgrTab = *(aPropagCoef.background());
	if(!fRcel)
		Init(backgrTab.Kmin(), backgrTab.Kmax());

	InitCoef(backgrTab);
	const CLinearFunc* backgr = ((CBackgroundTable&)backgrTab).getBackground();
	double rMinCEL = pow(10.,-ICS_CEL_SPLITTED_MAX_BINS_DIFF/BC().s);
	fa_discr.reset();
	fb_discr.reset();
	const int nn = Ranges().nE();
	const int mm = backgrTab.nK();
	int j = 0;
	const CBinning& midE = Ranges().midE();
	const CBinning& leftE = Ranges().E();
	const CBinning& midK = backgrTab.midBackgroundE();
	RKern kern(PICSeR, 0);
	for(j=0; j<nn; j++)
	{
		double E = midE[j];
		AKern discrAKern(fRdiscr, E);
		double aDiscr = backgr->integrate(discrAKern);
		ASSERT(aDiscr >= 0);
		fa_discr[j] = aDiscr;
		if(aDiscr > 0)
		{
			for(int i=0; i<j; i++)
			{
				double r1 = leftE[i]/E;
				double r2 = BC().ss1*r1;
				if(r2>rMinCEL)
					r2 = rMinCEL;
				if(r1<r2)
				{
					double sum = 0.;
					for(int k=0; k<mm; k++)
					{
						kern.SetW(midK[k]*E);
						sum += backgrTab.Fmid(k)*FuncUtils::integrateBode(r1,r2, kern);
					}
					fb_discr[i][j] = sum;
				}
			}
		}
	}

	int jMax = 2*Ranges().nE() + 1;

	if(!fEnergyLossRate.length())
		fEnergyLossRate.create(jMax);
	else
		fEnergyLossRate.reset();

	double E=Ranges().Emin();
	for(j=0; j<jMax; j++, E*=BC().ss2)
	{
		RhoKern celRhoKern(fRcel, fRho, E);
		double eLossCEL = backgr->integrate(celRhoKern);//mean Eph/E for CEL part
		ASSERT(eLossCEL>0.);
		fEnergyLossRate[j] = eLossCEL;
	}
	if(COEF_TEST_ON && redshift.z()==0.)
	{
		const CBinning& midE = Ranges().midE();
		CFilePtr file(FopenL(plt_local_dir + "ICSsplit","wt"));
		fprintf(file, "#E[eV]\tEloss_cel[Mpc^-1]\tEloss_discr[Mpc^-1]\n");
		for(j=0; j<nn; j++)
		{
			E=midE[j];

		//discret part
			double discrEloss = 0.;
			if(fa_discr[j]>0)
			{
				for(int i = 0; i<=j; i++)
				{
					discrEloss += fb_discr[i][j] * (E - midE[i])/E;
				}
			}
			//rate in Mpc^-1
			fprintf(file, "%lg\t%lg\t%lg\n", (E*units.Eunit*1e6), (units.Mpc_cm /units.Lunit*fEnergyLossRate[2*j+1]), (units.Mpc_cm /units.Lunit*discrEloss));
		}
	}
}

}//namespace couplings {


