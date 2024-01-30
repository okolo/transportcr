#include "Synchrotron.h"
#include "Ranges.h"
#include "Medium.h"
#include "TimeZ.h"
#include "Units.h"
#include "TableReader.h"


namespace couplings
{

//Synchrotron::Data::~Data()
//{
//	free(x);
//	free(y);
//};
//
//int Synchrotron::Data::InitData(const char *_fileName,double _lbackground,double _rbackground)
//{
//	FILE *ff;
//	int i,accurasy;
//	double xmem;
//	i=-1;
//	accurasy=-5;
//	if((ff=fopen(_fileName,"rt"))==NULL)
//	{
//		LOG_ERROR2("can't open file ", _fileName);
//		return 1;
//	};
//	if((fscanf(ff," ac=%d",&accurasy)<1)||accurasy<2)
//	{
//Data_error: LOG_ERROR4("error in file ", _fileName, " line ", i+2);
//		return 1;
//	};
//	ac=accurasy;
//	if((x=(double*)calloc(ac+1,sizeof(double)))==NULL)
//	{
//Data_memory: std::string str = "failed to load ";
//		MemoryExit((str + _fileName).c_str());
//		return 1;
//	};
//	if((y=(double*)calloc(ac+1,sizeof(double)))==NULL)
//		goto Data_memory;
//	i=0;
//	if(fscanf(ff," %lg %lg",&xmem,y)<2) goto Data_error;
//	x[0]=xmem;
//	ymin=(ymax=y[0]);
//	for(i=1;i<ac;i++)
//	{
//		if((fscanf(ff," %lg %lg",x+i,y+i)<2)||(x[i]<=xmem)) goto Data_error;
//		if(y[i]>ymax) ymax=y[i];
//		if(y[i]<ymin) ymin=y[i];
//	};
//	kx=ac/(x[ac-1]-x[0]);
//	k=(y[ac-1]-y[0])/(x[ac-1]-x[0]);
//	lbackground=_lbackground;
//	rbackground=_rbackground;
//	x[ac]=x[ac-1]+1/kx;
//	y[ac]=_rbackground;
//	return 0;
//};
//int Synchrotron::Data::InitDataL(char *_fileName,double _lbackground)
//{
//	if(InitData(_fileName,_lbackground)) return 1;
//	y[ac]=(rbackground=y[ac-1]);
//	return 0;
//};
//int Synchrotron::Data::InitDataR(char *_fileName,double _rbackground)
//{
//	if(InitData(_fileName,0,_rbackground)) return 1;
//	lbackground=y[0];
//	return 0;
//};
//int Synchrotron::Data::InitDataC(char *_fileName)
//{
//	if(InitData(_fileName)) return 1;
//	y[ac]=(rbackground=y[ac-1]);
//	lbackground=y[0];
//	return 0;
//};
//double Synchrotron::Data::GetData(double _x)
//{
//	int i;
//	double tem;
//	i=(int)floor((_x-x[0])*kx);
//	if(i<0) return lbackground;
//	if(i>=ac) return rbackground;
//	if(x[i]<_x)
//		for(;x[i+1]<=_x;i++);
//	else
//		for(;x[i]>_x;i--);
//	tem=y[i]+(_x-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
//	return tem;
//};


const double Synchrotron::B_const=0.000023667271986223;//4/9*alpha^2
//Synchrotron::Data Synchrotron::fGData;
//double Synchrotron::G(double x)
//	{
//		if((x<=0.0)||(x>=30.0))
//			return(0.0);
//		if(!fGData.IsInitialized() )
//		{
//			if(fGData.InitData(DATA_DIR"G.dat")!=0)
//				ThrowError("Failed to initialize SynchroSpectrum");
//		}
//		return fGData.GetData(x);
//	}

void Synchrotron::SetBackgrounds(const Medium& aPropagCoef)
	{
		Coupling::SetBackgrounds(aPropagCoef);

		double curB = aPropagCoef.MagnField();
		if(curB==fPrevB)
			return;//the field hasn't changed
		Syn.reset();
		Syn_e.reset();
		aeSyn.reset();
		if(curB<=0.)
			return;

		const int nn = Ranges().nE();
		int i=0;
		int j=0;
		double Ee,Eph;
		double B__1=curB/Me2;
		B__1 *= (B_const*B__1);

		Syn[0][0] = 0.;
        CVector norm_factor(nn);
        CVector fracEnergyConserved(nn);
		for(j=0 ; j<nn ; j++ )
		{
			Ee = Ranges().midE()[j];
            double synE = SynchroEnergy(Ee, curB);
            double intGtot = fIntG->f(Ee/synE); // maximal photon energy can not be higher than Ee
            // fraction of secondary gamma energy within (Emin, Emax) range
            double _fracEnergyConserved = (intGtot - fIntG->f(Ranges().midE()[0]/synE))/intGtot;
            fracEnergyConserved[j] = _fracEnergyConserved;
            double Efin = Ee*BC().ss_1;  // we always shift to just one energy bin since otherwise energy binning accuracy is not working
			aeSyn[j]=B__1*Ee*Ee/(Ee-Efin);
			if(j>0)
				Syn_e[j-1][j] = aeSyn[j];

            double sumPhotonE = 0.;// total energy going to photons
			for(i=0; i<=j; i++)
			{
				Eph = Ranges().midE()[i];
				Syn[i][j]= SynchroSpectrum(Ee, Eph, curB) * (BC().ss2 - BC().ss_2);
                sumPhotonE += (Syn[i][j] * Eph);
			};
            if(sumPhotonE > 0){
                double sumEprim = aeSyn[j] * (Ee-Efin); // total energy lost by electrons
                double mult = sumEprim * _fracEnergyConserved / sumPhotonE;
                for (i = 0; i <= j; i++) {
                    Syn[i][j] *= mult;
                }
                norm_factor[j] = mult;
            }
            else{
                ASSERT(_fracEnergyConserved<1e-10);
                norm_factor[j] = 1.;
            }
		}
		{
			ofstream out;
			if(COEF_TEST_ON && redshift.z()==0.)
			{
				string outputFile = plt_local_dir + "Synchrotron";
				out.open(outputFile.c_str());
				out << "E[eV]\tEc[eV]\tenergy loss length [Mpc]\t<ae> [Mpc^-1]\t<sumE_gamma>/<E_e_loss>\tfracEnergyConserved\tnorm_factor" << endl;
			}
			int jMin = 0;//nn/20;//must be adjusted
			int i, j;
			for(j=jMin;j<nn;j++)
			{
                Ee = Ranges().midE()[j];
                double rate=B__1*Ee; // - 1/E dE/dt

				double sumGammaE = 0;
				for(i=0;i<=j;i++)
				{
					ASSERT(Syn[i][j]>=0.);
                    sumGammaE += Ranges().midE()[i] * (Syn[i][j]);
				};
                sumGammaE /= Ee * (1. - BC().ss_1);

				if(out.is_open())
				{
					out << Ee / units.eV << "\t" << SynchroEnergy(Ee, curB) / units.eV << "\t" <<  1./rate/units.Mpc;
                    out << "\t" << aeSyn[j]*units.Mpc << "\t" << sumGammaE / aeSyn[j] << "\t";
                    out << fracEnergyConserved[j] << "\t" << norm_factor[j];
				}
				if(out.is_open()) out << endl;
			}
			if(out.is_open())
				out.close();
		}
		fPrevB = curB;
	}

	Synchrotron::Synchrotron():
		fPrevB(0.),
		fG(new CTableReader(DATA_DIR"G.dat", 2)),
        fIntG((CLinearFunc*)fG.integral())
	{
		const int nn = Ranges().nE();
		Syn.create(nn);
		Syn_e.create(nn);
		aeSyn.create(nn);

		AddChannel(new Channel_e_e(this, EElectron));
		AddChannel(new Channel_e_e(this, EPositron));
		AddChannel(new Channel_e_gamma(this, EElectron));
		AddChannel(new Channel_e_gamma(this, EPositron));
	}

	void Synchrotron::Channel_e_e::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		int iPrim, iSec;
		for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		{
			aCoef.Add(iPrim, iPrim, -fCoupling.aeSyn[iPrim]);
			for(iSec = 0; iSec<=iPrim; iSec++)
				aCoef.Add(iSec, iPrim, fCoupling.Syn_e[iSec][iPrim]);
		}
	}

	void Synchrotron::Channel_e_gamma::Coef(CMatrixAddOnlyView& aCoef) const
	{
		const int nn = Ranges().nE();
		int iPrim, iSec;
		for(iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
		{
			for(iSec = 0; iSec<=iPrim; iSec++)
				aCoef.Add(iSec, iPrim, fCoupling.Syn[iSec][iPrim]);
		}
	}
}
