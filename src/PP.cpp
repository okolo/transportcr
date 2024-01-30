#include "PP.h"

#include "Ranges.h"
#include "Medium.h"
#include "SafeOutput.h"
#include "TimeZ.h"
#include "Units.h"

namespace couplings
{
	PP::PP(void):
		cs(ParticleData::getParticleMass(EElectron))
	{
		AddChannel(new Channel_gamma_e(this, EElectron));
		AddChannel(new Channel_gamma_e(this, JOIN_ELECTRON_POSITRON ? EElectron : EPositron));
		AddChannel(new Channel_gamma_gamma(this));
	}

	void PP::SetBackgrounds(const Medium& aPropagCoef)
	{
		const bool LEE_ADJUSTMENT_PP = false;//coef adjustments (experimental)
		Coupling::SetBackgrounds(aPropagCoef);

		const int nn = Ranges().nE();
		if(!aph.length())
		{
			ce.create(nn);
			aph.create(nn);
			intPrdrPPlow.create(nn);
			intPdrPPlow.create(nn);
		}
		else
		{
			ce.reset();
			aph.reset();
			intPrdrPPlow.reset();
			intPdrPPlow.reset();
		}
		InitCoef(*(aPropagCoef.background()));
		if(LEE_ADJUSTMENT_PP)
			AdjustCoef();
		FinalCoefAdjustment();
	}

	void PP::FinalCoefAdjustment()
	{
		const int nn = Ranges().nE();
		for(int j=0;j<nn;j++)
		{
			double sum2=0;
			for(int i=0;i<=j;i++)
			{
				sum2+=ce[i][j];
			}
			if(sum2>aph[j])
			{
				if((sum2-aph[j])/sum2>1e-10)
				{
					ASSERT(false);
					LOG_ERROR_ONCE("aph[j]<sum(ce[i][j])");
				}
				aph[j]=sum2;
			}
		}
	}

	void PP::Channel_gamma_gamma::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		int iPrim;
		for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
			aCoef.Add(iPrim, iPrim, -fCoupling.aph[iPrim]);
	}

	void PP::DebugOutput(const char* aFolder) const
	{
		ASSERT(aph.length()>0);//coefficients must be initialized
		const int nn = Ranges().nE();
		double ycoeff=units.Lunit/units.Mpc_cm;
		double E = Ranges().Emin()*units.Eunit/units.outEunitMeV;
		CFilePtr eInt(Fopen("ppIntLengthC", aFolder, "wt"));
		for(int i=0; i<nn; i++,E*=BC().ss1)
		{
			if(aph[i]>0)
					fprintf(eInt,"%lg\t%lg\n", E, ycoeff/aph[i]);
		}
	}

	void PP::Channel_gamma_e::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		int iPrim, iSec;
		for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		{
			for(iSec = 0; iSec<=iPrim; iSec++)
				aCoef.Add(iSec, iPrim, fCoupling.ce[iSec][iPrim]);
		}
	}
	void PPapr10::InitCoef(const CBackgroundTable& aBackgr)
	{
		const int nn = Ranges().nE();
		const int mm = aBackgr.nK();
		const CBinning& midE = Ranges().midE();

		const double Emin = Ranges().Emin();

		for(int j=0; j<nn; j++)
		{
			double E = midE[j];
			for(int k=0; k<mm; k++)
			{
				double Fb=aBackgr.Fmid(k);
				double b=aBackgr.midBackgroundE()[k];

				double r = Emin/E;
				for(int i=0; i<=j; i++, r*=BC().ss1)
				{
					ce[i][j]+=Fb*cs.PPPdrIntegral(E, b, r, r*BC().ss1);
				}
			}
		}
		for(int k=0; k<mm; k++)
		{
			double Fb=aBackgr.Fmid(k);
			double b=aBackgr.midBackgroundE()[k];
			for(int i=0; i<nn; i++)
			{
				double E = midE[i];
				aph[i]+=Fb*cs.RPP(E,b);
				intPrdrPPlow[i] += Fb*cs.lowerPPPrdrIntegral(E,b,Emin/E);
				intPdrPPlow[i] += Fb*cs.lowerPPPdrIntegral(E,b,Emin/E);
			}
		}
	}

	void PPapr10::AdjustCoef()
	{
		//double usageLimit = 0.015;//maximal descrepancy between total sigma and integral of dif sigma for energy conservation based adjustment of coefficients
		CSafeOutput ppOut;
		if(COEF_TEST_ON)
		{
			ppOut.open(plt_local_dir + "PPinfoZ" + ToString(redshift.z()));
			ppOut << "E\tEsec\tratio\t<adjustment status>\n" ;
		}

		const int nn = Ranges().nE();
		int iMin = 1;
		const double binDifRatio = 1.-BC().ss_1;
		
		int i,j;
		
		//PP
		bool secondAttempt = false;
		for(i=iMin; i<nn; i++)
		{
			double Ei = Ranges().midE()[i];
			ppOut << i << "\t" << Ei*units.Eunit*1e6 << "\t";
			ASSERT(aph[i]>=0);
			
			if(aph[i]==0.)
			{
				ppOut << "no interaction" << endl;
				continue;
			}
			
			double sumEP = intPrdrPPlow[i];
			double sumP = intPdrPPlow[i];
			int k, l;
			for(k=i; k>=0 && ce[k][i]==0; k--);//finding highest secondary energy bin
			for(l=0; l<=i && ce[l][i]==0.; l++);//finding lowest secondary energy bin
			double Ek = Ranges().midE()[k];
			//bool always = true;
			if(sumP>0)
			{//adjust ce[k-1][i] and ce[k][i]
				ASSERT(l==0);
				//double Ei_1 = Ranges().midE()[i-1];
				int jMax = k-1;
				
				for(j=0; j<jMax; j++)		
				{
					if(ce[j][i]>0)
					{
						sumEP += Ranges().midE()[j]*ce[j][i];
						sumP += ce[j][i];
					}
					else
					{
						ASSERT(ce[j][i]==0.);
					}
				};
				
				double ce_k_1 = (aph[i]*(Ek-0.5*Ei) - Ek*sumP + sumEP)/(Ek*binDifRatio);
				double ce_k = aph[i] - sumP - ce_k_1;
				
				if(ce_k_1>=0&&ce_k>=0)
				{
					secondAttempt = false;
					ppOut << "PP upper coef adjustment OK";
					ce[k-1][i] = ce_k_1;
					ce[k][i] = ce_k;
				}
				else if(ce_k<0)
				{
					if(secondAttempt)
					{
						const char* errMsg = "PP upper coef adjustment fatal error: second attempt failed";
						ppOut << errMsg;
						ASSERT(false);
						throw errMsg;
					}
					ppOut << "PP upper coef adjustment second attempt\n";
					ce[k][i] = 0.;
					i--;
					secondAttempt = true;
					continue;
					//old adjustment was
					//aph[i] = sumP + ce[k-1][i];
					//ppOut << "new mode FAILED";
				}else
				{
					const char* errMsg = "PP upper coef adjustment fatal error: ce_k_1<0";
					ppOut << errMsg;
					ASSERT(false);
					throw errMsg;
				}
			}
			else
			{//adjust ce[l][i] and ce[k][i]
				for(j=l+1; j<k; j++)		
				{
					ASSERT(ce[j][i]>0.);
					sumEP += Ranges().midE()[j]*ce[j][i];
					sumP += ce[j][i];
				};
				
				double ce_l = (aph[i]*(2.*Ek-Ei) - 2.*Ek*sumP + 2.*sumEP)/(2.*Ek - Ranges().midE()[l]);
				double ce_k = aph[i] - sumP - ce_l;
				
				if(ce_l>=0&&ce_k>=0)
				{
					secondAttempt = false;
					ppOut << "PP lower coef adjustment OK";
					ce[l][i] = ce_l;
					ce[k][i] = ce_k;
				}
				else
				{
					if(secondAttempt)
					{
						const char* errMsg = "PP lower coef adjustment fatal error: second attempt failed";
						ppOut << errMsg;
						ASSERT(false);
						throw errMsg;
					}
					ppOut << "PP lower coef adjustment second attempt\n";
					ce[k][i] = 0.;
					i--;
					secondAttempt = true;
					continue;
				}			
			}
			ppOut << "\ti-k=" << i-k << "\tlower_a/a="<< intPdrPPlow[i]/aph[i] << "\tlowest secondary bin "<< l <<  "\n";
		};
		ppOut.close();
	}

	void PPold::InitCoef(const CBackgroundTable& aBackgr)
	{
		int i,j,k;
		double Emin = Ranges().Emin();
		const int nn = Ranges().nE();
		const int mm = aBackgr.nK();
		//const double bmin = aBackgr.Kmin();

		double tem1, tem2, E, E1, sss;
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
				const double Fb = aBackgr.Fmid(k);
				const double b = aBackgr.midBackgroundE()[k];

				tem2=cs.PPPe(E,b,Emin);
				VERIFY(tem2>=0);
				for(E1=Emin*BC().ss1,i=0;i<j;i++,E1*=BC().ss1)
				{
					tem1=cs.PPPe(E,b,E1);
					ce[i][j]+=E1*sss*(a1*tem1+a1*tem2+a2*(cs.PPPe(E,b,x2*E1)+cs.PPPe(E,b,x4*E1))+a3*cs.PPPe(E,b,x3*E1))*Fb;
					tem2=tem1;
				};
				tem1=cs.PPPe(E,b,0.99*E);
				ce[j][j]+=E*(1.0-BC().ss_2)*(a1*tem1+a1*tem2+a2*(cs.PPPe(E,b,xx2*E)+cs.PPPe(E,b,xx4*E))+a3*cs.PPPe(E,b,xx3*E))*Fb;
			}
		}
		for(k=0 ;k<mm; k++)
		{
			const double Fb = aBackgr.Fmid(k);
			const double b = aBackgr.midBackgroundE()[k];
			for(i=0,E=Emin*BC().ss2;i<nn;i++,E*=BC().ss1)
			{		
				aph[i]+=Fb*cs.RPP(E,b);
				intPrdrPPlow[i] += Fb*cs.lowerPPPrdrIntegral(E,b,Emin/E);
				intPdrPPlow[i] += Fb*cs.lowerPPPdrIntegral(E,b,Emin/E);
			}
		}
	}

	void PPold::AdjustCoef()
	{
		//double usageLimit = 0.015;//maximal descrepancy between total sigma and integral of dif sigma for energy conservation based adjustment of coefficients
		CSafeOutput ppOut;
		if(COEF_TEST_ON)
		{
			ppOut.open(plt_local_dir + "PPinfoZ" + ToString(redshift.z()));
			ppOut << "E\tEsec\tratio\t<adjustment status>\n" ;
		}

		const int nn = Ranges().nE();
		int iMin = 1;
		const double binDifRatio = 1.-BC().ss_1;
		
		int i,j;
		
		for(i=iMin; i<nn; i++)
		{
			double Ei = Ranges().midE()[i];
			ASSERT(aph[i]>=0);
			if(aph[i]==0.)
			{
				ppOut << i << "\t" << Ei*units.Eunit*1e6 << "\t no interaction" << endl;
				continue;
			}
			//double Ei_1 = Ranges().midE()[i-1];
			double sumEP = intPrdrPPlow[i];
			double sumP = intPdrPPlow[i];
			int k;
			for(k=i; k>=0 && ce[k][i]==0; k--);//finding highest secondary energy bin
			int jMax = k-1;
			double Ek = Ranges().midE()[k];
			int l = 0;//lowest secondary bin
			for(j=0; j<jMax; j++)		
			{
				if(ce[j][i]>0)
				{
					sumEP += Ranges().midE()[j]*ce[j][i];
					sumP += ce[j][i];
				}
				else
				{
					l++;
					ASSERT(ce[j][i]==0.);
				}
			};
			ppOut << i << "\t" << Ei*units.Eunit*1e6 << "\t";
		

			double ce_k_1 = (aph[i]*(Ek-0.5*Ei) - Ek*sumP + sumEP)/(Ek*binDifRatio);
			double ce_k = aph[i] - sumP - ce_k_1;
			
			if(ce_k_1>=0&&ce_k>=0)
			{
				ppOut << "PPold mode OK";
				ce[k-1][i] = ce_k_1;
				ce[k][i] = ce_k;
			}
			else if(ce_k<0)
			{
				ce[k][i] = 0.;
				aph[i] = sumP + ce[k-1][i];
				ppOut << "PPold mode FAILED";
			}else
			{
				ppOut << "PPold mode FATAL ERROR!";
				ASSERT(false);
			}

			ppOut << "\ti-k=" << i-k << "\tlower_a/a="<< intPdrPPlow[i]/aph[i] << "\tlowest secondary bin "<< l <<  "\n";
		};
		ppOut.close();
	}
}

