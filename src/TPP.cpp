#include "TPP.h"
#include "Ranges.h"
#include "Medium.h"
#include "Units.h"
#include "TimeZ.h"

namespace couplings
{
/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
/* Triplet Pair Production angle averaged full cross section */
	double TPP::RTPP(double E,double b)
	{

		double s,y;
		s=4*E*b;
		if(s<13.4*Me2)
			return(0.0);
		y=((151*Me2-65*s)/3+7*(s-2*Me2)*log(s/Me2)+140.842*Me2*Me2/s)/s;
		y*=1.72708440115e-7/Me2;
		return(y);
	};

/*88888888888888888888888888888888888888888888888888888888888888888888888888888888888
88888           Continuous energy loss TPP approximation           888888888888888888
8888888888888     (-1/E)dE/dt=Integral(n(b)F_TPP(E,b))db           88888888888888888*/
	double TPP::F_TPP(double E,double b)
	{
		const double A=1.908427e-08;//1.768/36*alpha^3

		const double a1=369.33333;// 1108/3
		const double a2=-38.026667;// -2852/75
		const double a3=11.2;//      56/5

		const double C0=198.860592333824;//a1*x^(1/4)+a2*x^(5/4)+a3*(x-5.0)*x^(1/4)*log(x) , where x=13.4

		double x=1+4*E*b/Me2;// S_max/Me^2

		if(x<13.4)// from the condition sigma>0
			return 0.0;

		double x4=pow(x,0.25);

		double result=a1*x4+a2*x*x4+a3*(x-5.0)*x4*log(x)-C0;
		result*=A*Me2/b/E/b/E;
		return result;
	};

	/*8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888*/
	  /*  Triplet Pair Production angle averaged differential cross section
	   (per e energy) */
	double TPP::PTPP(double E,double b,double nE)
	/* E-cosmic ray electron energy
	   b-background photon energy
	   nE-the energy of one of the particles of the produced pair */
	{
		double Emax,Emin,s,x;
		s=E*b+Me2;
		if(s<13.4*Me2)
			return(0.0);
		Emax=(Emin=1.0-0.5*Me2/s);
		x=sqrt(1.0-2.0*Me2/s);
		Emax+=x;
		Emin-=x;
		x=2.0*E*s/(4.0*s+Me2);
		Emax*=x;
		Emin*=x;
		if((nE>Emax)||(nE<Emin))
			return(0.0);
		//if(TPP_ENERGY_CONSERVATION)
			x=0.125*F_TPP(E,b)*E/(pow(Emax,0.25)-pow(Emin,0.25));
		//else
		//	x=RTPP(E,b)*0.75/(pow(Emin,-0.75)-pow(Emax,-0.75));
		return x*pow(nE,-1.75);
	};

	TPP::TPP(void)
	{
		AddChannel(new Channel_e_e(this, EElectron));
		AddChannel(new Channel_e_e(this, EPositron));
		AddChannel(new Channel_e_ae(this, EElectron, JOIN_ELECTRON_POSITRON ? EElectron : EPositron));
		AddChannel(new Channel_e_ae(this, EPositron, EElectron));
	}

	void TPP::SetBackgrounds(const Medium& aPropagCoef)
	{
		Coupling::SetBackgrounds(aPropagCoef);
		const SmartPtr<CBackgroundTable>& ebl = fBackground->background();
		const int nn = Ranges().nE();
		const int mm = aPropagCoef.background()->nK();
		const CBinning& midK = aPropagCoef.background()->midBackgroundE();
		const CBinning& leftE = Ranges().E();
		const CBinning& midE = Ranges().midE();
		const double ssa=0.5*(1.0+BC().ss_1);
		const double sss6=(1.0-BC().ss_1)/6.0;
		if(f_a.length()==0)
		{
			f_a.create(nn);
			f_b.create(nn);
		}
		else
		{
			f_a.reset();
			f_b.reset();
		}
		for(int k=0; k<mm; k++)
		{
			double b = midK[k];
			double Fb=ebl->Fmid(k);
			for(int i=0; i<nn; i++)
			{
				f_a[i]+=Fb*F_TPP(leftE[i],b);
				ASSERT(f_a[i]>=0);
			};
		};

		for(int j=0; j<nn; j++)
		{
			double E = midE[j];
			for(int k=0; k<mm; k++)
			{
				double b = midK[k];
				double Fb=ebl->Fmid(k);
				double tem2=PTPP(E,b,Ranges().Emin());
				ASSERT(tem2>=0);
				for(int i=0; i<j; i++)
				{
					double E1 = leftE[i+1];
					double tem1=PTPP(E,b,E1);
					f_b[i][j] += sss6*E1*(tem1+tem2+4.0*PTPP(E,b,ssa*E1))*Fb;
					tem2=tem1;
				};
			};
		};
	}

	void TPP::Channel_e_e::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		double CELconst = 1/(BC().ss1-1.0);

		for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		{
			//CEL
			aCoef.Add(iPrim, iPrim, -fCoupling.f_a[iPrim]*CELconst);
			if(iPrim<nn-1)
				aCoef.Add(iPrim, iPrim+1, fCoupling.f_a[iPrim+1]*CELconst);

			for(int iSec = 0; iSec<=iPrim; iSec++)
				aCoef.Add(iSec, iPrim, fCoupling.f_b[iSec][iPrim]);
		}
	}

	void TPP::Channel_e_ae::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
			for(int iSec = 0; iSec<=iPrim; iSec++)
				aCoef.Add(iSec, iPrim, fCoupling.f_b[iSec][iPrim]);
	}

	void TPP::DebugOutput(const char* aFolder) const{
		const int nn = Ranges().nE();
		std::string path = aFolder;
		path += "TPP";

		if(redshift.z()>0)
			path = path + "_z" + ToString(redshift.z());

		ofstream file(path.c_str());
		const CBinning& E=Ranges().E();//using Ranges().E and not Ranges().midE for CEL approximation
		const double p_mass = ParticleData::getParticleMass(EProton);
		file << "# E [eV]   1/E*dE/dt [Mpc^-1]\n";

		for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		{
			if(f_a[iPrim]>0){
				file << (E[iPrim]/units.eV) << "\t" << f_a[iPrim]*units.Mpc << "\n";
			}
		}
		file.close();
	}
}


