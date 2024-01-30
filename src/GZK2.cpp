#include <math.h>
#include <fstream>
#include "GZK2.h"
#include "Ranges.h"
#include "Medium.h"
#include "Nucleus.h"
#include "TimeZ.h"
#include "TableReader.h"
#include "Units.h"
#include <limits>

using namespace std;

namespace couplings{

const TParticle GZK::secondaries[] =
{
	EElectron,
	EPositron,
	EPhoton,
	ENeutrinoE,
	ENeutrinoM,
	ENeutrinoAE,
	ENeutrinoAM,
	ENeutron,
	EProton,
	EEndParticle
};

const TParticle GZK::primaries[] =
{
	ENeutron,
	EProton,
	EEndParticle
};
bool GZK::NO_PI_PRODUCTS = false;//ignore input of photopions to leptons and photons spectra
bool GZK::NO_PI_NUCLEI = false;//do not take into account pion production for nuclei

GZK::GZK()
{
	READ_BOOL_SETTING(NO_PI_PRODUCTS);
	READ_BOOL_SETTING(NO_PI_NUCLEI);

	for(int prim = 0; prim<=1; prim++)
	{
		fSecDistributions[prim] = 0;
		fSigmas[prim] = 0;
	}

	const CParticleList& pf = *CParticleList::Instance();
	if(JOIN_ELECTRON_POSITRON)
	{
		throw "JOIN_ELECTRON_POSITRON mode is not supported yet by GZK coupling";
	}
	AddChannel(new Channel_N_x(this, EProton, EProton));
	AddChannel(new Channel_N_x(this, EProton, ENeutron));
	AddChannel(new Channel_N_x(this, ENeutron, ENeutron));
	AddChannel(new Channel_N_x(this, ENeutron, EProton));
	if(!NO_PI_PRODUCTS)
	{
		for(int i = 0; secondaries[i]<ENeutron; i++)
		{
			TParticle p = secondaries[i];
			if(!pf.IsEnabled(p))
				continue;
			AddChannel(new Channel_N_x(this, EProton, p));
			AddChannel(new Channel_N_x(this, ENeutron, p));	
		}
	}
	if(NO_PI_NUCLEI)
		return;
	FOR_ALL_NUCLEI_INVOLVED(particle)
	{
		AddChannel(new Channel_A_A(this, particle));
		
		if(!NO_PI_PRODUCTS)
		{
			for(int i = 0; secondaries[i]<ENeutron; i++)
			{
				TParticle p = secondaries[i];
				if(pf.IsEnabled(p))
					AddChannel(new Channel_A_x(this, particle, p));
			}
		}
	}
}

GZK::~GZK()
{
}

void GZK::DebugOutput(const char* aFolder) const
{
//print energy loss rate -dE/dt for protons due to GZK assuming instant beta decay of neutrons
//and neglecting energy loss in beta-decay
	const int nn = Ranges().nE();
	string path = aFolder;
	path = path + "GZK";
	Mkdir(path);
	path = path + DIR_DELIMITER_STR + ToString(redshift.z());
	ofstream file(path.c_str());
	const CBinning& E=Ranges().midE();
	int index = ((std::map<int, int>&)fMapP)[EProton];
	const CMatrix& pp = fP[1](index);
	index = ((std::map<int, int>&)fMapP)[ENeutron];
	const CMatrix& pn = fP[1](index);
	double dt_dz = CTimeZ::diffT(redshift.z());

	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		double a = fR[1][iPrim];
		double eLoss=0;
		double eLossZ=0;
		if(a!=0)
		{
			double sum=0.;
			for(int iSec = 0; iSec <= iPrim; iSec++)
				sum +=  E[iSec]*(pp[iSec][iPrim] + pn[iSec][iPrim]);
			eLoss = (a*E[iPrim]-sum)*units.Eunit*1e6/units.Lunit*units.Mpc_cm;//[eV/Mpc]
			eLossZ = eLoss*units.Lunit/units.Mpc_cm *dt_dz;// dE/dz [eV]
		}
		file << (E[iPrim]*units.Eunit*1e6) << "\t" << eLossZ << "\t" << eLoss << "\n";
	}
	file.close();
}

void GZK::init(double aMaxBackgroundE)
{
	const CParticleList& pf = *CParticleList::Instance();
	const int nn = Ranges().nE();
	std::string dataDir = DATA_DIR "sophia2" DIR_DELIMITER_STR;
	std::string secDataDir = dataDir + ToString(BC().NoOfBinsPerDecade);

	if(!directory_exists(secDataDir))
		ThrowError(secDataDir + " directory not found\nGenerate SOPHIA tables for " + ToString(BC().NoOfBinsPerDecade) + " bins per decade or change binning!");

	secDataDir = secDataDir + DIR_DELIMITER_STR;

	double kMaxPossible = Ranges().Emax()/ParticleData::getParticleMass(EProton)*2.*aMaxBackgroundE;
	for(int prim=0; prim<=1; prim++)
	{
		const char* primName = ParticleData::getParticleFileName(primaries[prim]);
		CTableReader* reader = new CTableReader(dataDir + primName, 2);
		CVarVector& k = reader->getColumn(0);//in GeV
		CVarVector& sigma = reader->getColumn(1);//in mubarn
		for(int j=k.length()-1; j>=0; j--)
		{
			//converting to internal units
			k[j] *= (1e3/units.Eunit);
			sigma[j] *= (1e-6*units.barn);
		}

		fSigmas[prim] = new CLogScaleLinearFunc(reader);
		fSigmas[prim]->SetExtension(ExtConst, false, true);
		const double csTableKmax = k[k.length()-1];
		if(csTableKmax<kMaxPossible)
		{
			ASSERT(false);
			std::string message = "SOPHIA total cross section tables have Kmax/GeV=" +
					ToString(csTableKmax*1e-3*units.Eunit) +	" which is less than required value " +
					ToString(kMaxPossible*1e-3*units.Eunit) + " -> extrapolation will be used";
			WARN_ONCE(message);
		}

		for(int sec = 0; secondaries[sec]!=EEndParticle; sec++)
		{
			TParticle sp = secondaries[sec];
			if(sp<ENeutron && (NO_PI_PRODUCTS || !pf.IsEnabled(sp)))
				continue;
			fP[prim].add(new CMatrix(nn));
			std::string secdFile = secDataDir + primName + "-" + ParticleData::getParticleFileName(sp);
			SecDistribution* secd = new SecDistribution(secdFile);
			fSecDistributions[prim].add(secd);
			if(secd->Primary()!=primaries[prim])
				ThrowError("Unexpected primary type in file " + secdFile);
			if(secd->Secondary()!=sp)
				ThrowError("Unexpected secondary type in file " + secdFile);
			if(secd->NbinsPerDecade()!=BC().NoOfBinsPerDecade)
				ThrowError("Unexpected NbinsPerDecade=" + ToString(secd->NbinsPerDecade()) + " in file " + secdFile);
			if(secd->Kmax()<kMaxPossible)
			{
				ASSERT(false);
				std::string message = "SOPHIA diff cross section tables have Kmax/GeV=" +
						ToString(secd->Kmax()*1e-3*units.Eunit) +	" which is less than required value " +
						ToString(kMaxPossible*1e-3*units.Eunit) + " -> extrapolation will be used";
				WARN_ONCE(message);
			}

			if(prim==0)
				fMapP[sp] = fP[prim].length()-1;
		}
		fR[prim].create(nn);
		fCEL[prim].create(nn);
	}
}

void GZK::SetBackgrounds(const Medium& aPropagCoef)
{
	const CParticleList& pf = *CParticleList::Instance();
	const int nn = Ranges().nE();

	if(fSigmas[0].isNull())//lazy init
	{
		init(aPropagCoef.background()->Kmax());
	}
	else
	{
		for(int prim=0; prim<=1; prim++)
		{
			for(int i=fP[prim].length()-1; i>=0; i--)
				fP[prim](i).reset();
			fR[prim].reset();
			fCEL[prim].reset();
		}
	}
	const CBackgroundIntegral& backgrI = aPropagCoef.BackgroundIntegral();
	for(int prim=0; prim<=1; prim++)
	{
		CLogScaleLinearFunc* totSigma = fSigmas[prim];
		double kThr = totSigma->firstX();//threshold photon energy in nucleon rest frame
		double m = ParticleData::getParticleMass(primaries[prim]);
		double Eth = 0.5*kThr/backgrI.Kmax()*m;//minimal nucleon lab frame energy
		if(Eth>=Ranges().Emax())
			continue;//no interactions (below threshold)
		int iMin = Ranges().E().findLeftX(Eth);
		if(iMin<0)//Eth<Emin
			iMin = 0;
		for(int iEprim=iMin; iEprim<nn; iEprim++)
		{
			double gamma = Ranges().midE()[iEprim]/m;
			double maxK = 2.*backgrI.Kmax()*gamma;
			fR[prim][iEprim] = backgrI.calculateRlogsc(gamma, totSigma, kThr, maxK, 1e-3);

			for(int i = 0; secondaries[i]<EEndParticle; i++)
			{
				TParticle sec = secondaries[i];
				if(sec<ENeutron && (NO_PI_PRODUCTS || !pf.IsEnabled(sec)))
					continue;

				int secIndex = fMapP[sec];
				CMatrix& P = fP[prim](secIndex);
				SecDistribution& secDistr = fSecDistributions[prim](secIndex);
				//SecKern secKern(secDistr, *totSigma);
				int iEsecMin = iEprim - secDistr.MaxDeltaE();
				int iEsecMax = iEprim - secDistr.MinDeltaE();
				if(iEsecMin>iEsecMax)//zero diff cross section
				{
					ASSERT(false);//check table reading procedure or get rid of unnecessary tables
					continue;
				}
				for(int iEsec = iEsecMin<0?0:iEsecMin; iEsec<=iEsecMax; iEsec++)
				{
					//secKern.SetDeltaE(iEprim-iEsec);
					double R = fR[prim][iEprim];
					ASSERT(R>=0);
					if(R>0)
					{
						//double abserr = 1e-6*R;
						//double val1 = backgrI.calculateRlogsc(gamma, &secKern, kThr, maxK, 1e-3, abserr);
						double val2 = secDistr.CalculateRate(backgrI, *totSigma, iEprim-iEsec, gamma);
						/*
						if(val1>0 || val2>0)
						{
							double err = fabs(val1-val2)/(val1+val2)*2.;
							ASSERT(err<1e-2);
						}*/
						P[iEsec][iEprim] = val2;
					}
				}
			}
		}
	}

	if(!NO_PI_NUCLEI) 
	{
		double ss1 = BC().ss1;
		double coef_CEL = ss1/(ss1-1.);
		CAutoDeletePtrArray<CMatrix>& fPn = fP[0];
		CAutoDeletePtrArray<CMatrix>& fPp = fP[1];
		int indexP = fMapP[EProton];
		int indexN = fMapP[ENeutron];
		for(int j=0; j<nn; j++)
		{
			double E=Ranges().midE()[j];
			double ann = 0.;
			double app = 0.;
			double anp = 0.;
			double apn = 0.;
			double bnn = 0.;
			double bpp = 0.;
			double bnp = 0.;
			double bpn = 0.;
			for(int i=0; i<j; i++)
			{
				double E1=Ranges().midE()[i];
				bnn += E1*fPn(indexN)[i][j];
				bpp += E1*fPp(indexP)[i][j];
				bnp += E1*fPn(indexP)[i][j];
				bpn += E1*fPp(indexN)[i][j];
				ann += fPn(indexN)[i][j];
				app += fPp(indexP)[i][j];
				anp += fPn(indexP)[i][j];
				apn += fPp(indexN)[i][j];
			}
			double aPhotopion_nn = (ann - bnn/E)*coef_CEL;
			double aPhotopion_pp = (app - bpp/E)*coef_CEL;
			double aPhotopion_np = (anp - bnp/E)*coef_CEL;
			double aPhotopion_pn = (apn - bpn/E)*coef_CEL;
			fCEL[0][j] = aPhotopion_nn + aPhotopion_np;
			fCEL[1][j] = aPhotopion_pp + aPhotopion_pn;
		}
	}
}

void GZK::GetMultGZKcel(double& multN, double& multP, TParticle aPrim)
{
	double Z = CNucleus::getZ(aPrim);
	double A = CNucleus::getA(aPrim);

	double photopionMultiplierRN = pow((double)A,-1./3.);//	= A^{2/3}(cross section)  / A(energy loss suppression)
	multN = ((double)(A-Z))*photopionMultiplierRN/((double)A);//  * N/A (fraction of neutrons)
	multP = ((double)Z)*photopionMultiplierRN/((double)A);//  * Z/A (fraction of protons)
}

void GZK::GetMultGZKsec(double& multN, double& multP, TParticle aPrim)
{
	double Z = CNucleus::getZ(aPrim);
	double A = CNucleus::getA(aPrim);

	double photopionMultiplierRN = pow((double)A,2./3.);//	= A^{2/3}(cross section)
	multN = ((double)(A-Z))*photopionMultiplierRN/((double)A);//  * N/A (fraction of neutrons)
	multP = ((double)Z)*photopionMultiplierRN/((double)A);//  * Z/A (fraction of protons)
}

void GZK::Channel_N_x::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	int iPrim;
	if(Primary()==Secondary())
	{
		const CVector& R =  Primary()==EProton ? fCoupling.fR[1] : fCoupling.fR[0];
		for(iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
			aCoef.Add(iPrim, iPrim, -R[iPrim]);
	}
	const int index = fCoupling.fMapP[Secondary()];
	const CMatrix& P = Primary()==EProton ? fCoupling.fP[1](index) : fCoupling.fP[0](index);
	for(iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
		for(int iSec = 0; iSec <= iPrim; iSec++)
			aCoef.Add(iSec, iPrim, P[iSec][iPrim]);
}

void GZK::Channel_A_A::Coef(CMatrixAddOnlyView& aCoef) const
{
	double multN, multP;
	GZK::GetMultGZKcel(multN, multP, Primary());
	const int nn = Ranges().nE();

	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, - multN * fCoupling.fCEL[0][iPrim]);
		aCoef.Add(iPrim, iPrim, - multP * fCoupling.fCEL[1][iPrim]);
		if(iPrim<nn-1)
		{
			aCoef.Add(iPrim, iPrim+1, multN * fCoupling.fCEL[0][iPrim+1]);
			aCoef.Add(iPrim, iPrim+1, multP * fCoupling.fCEL[1][iPrim+1]);
		}
	}
}

void GZK::Channel_A_x::Coef(CMatrixAddOnlyView& aCoef) const
{
	double multN, multP;
	GZK::GetMultGZKsec(multN, multP, Primary());

	const int nn = Ranges().nE();

	int iPrim;
	const int index = fCoupling.fMapP[Secondary()];

	const CMatrix& Pp = fCoupling.fP[1](index);
	const CMatrix& Pn = fCoupling.fP[0](index);

	for(iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
		for(int iSec = 0; iSec <= iPrim; iSec++)
		{
			aCoef.Add(iSec, iPrim, multP*Pp[iSec][iPrim]);
			aCoef.Add(iSec, iPrim, multN*Pn[iSec][iPrim]);
		}
}


GZK::SecDistribution::SecDistribution(const std::string& aPath, bool aBinaryMode):
fPrimary(EEndParticle),
fSecondary(EEndParticle),
fKmin(-1.),
fKmax(-1),
fStep(-1),
fNbinsPerDecade(-1),
fNtrials(-1),
fGlobalMinDeltaE(std::numeric_limits<int>::max()),
fGlobalMaxDeltaE(std::numeric_limits<int>::min())
{
	ios_base::openmode openMode = aBinaryMode?ios::in|ios::binary : ios::in;
	ifstream inFile(aPath.c_str(), openMode);

	int part = -1;
	readL(inFile, part, aBinaryMode);
	ASSERT(part>=0 && part<EEndParticle);
	fPrimary = (TParticle)part;
	readL(inFile, part, aBinaryMode);
	ASSERT(part>=0 && part<EEndParticle);
	fSecondary = (TParticle)part;

	readL(inFile, fKmin, aBinaryMode);
	readL(inFile, fNbinsPerDecade, aBinaryMode);
	fStep = pow(10.,1./fNbinsPerDecade);
	readL(inFile, fNtrials, aBinaryMode);
	double norm = 1./fNtrials;
	fKmin *= (1e3/units.Eunit);
	int minDeltaE,maxDeltaE;
	while(read(inFile, minDeltaE, aBinaryMode))
	{
		readL(inFile, maxDeltaE, aBinaryMode);
		std::vector<double>* secondaries = new std::vector<double>();
		fData.add(secondaries);
		if(minDeltaE>=0)//secondaries exist
		{
			ASSERT(maxDeltaE>=0);
			if(fGlobalMinDeltaE>minDeltaE)
				fGlobalMinDeltaE = minDeltaE;
			if(fGlobalMaxDeltaE<maxDeltaE)
				fGlobalMaxDeltaE = maxDeltaE;
			unsigned long nEvents = 0;
			for(int i=minDeltaE; i<=maxDeltaE; i++)
			{
				readL(inFile, nEvents, aBinaryMode);
				ASSERT(nEvents>=0);
				secondaries->push_back(norm*nEvents);
			}
		}
		fMinDeltaE.push_back(minDeltaE);
	}
	if(!inFile.eof())
	{
		ThrowError("Data format error in " + aPath);
	}
	fKmax = fKmin*pow(10.,((double)(fMinDeltaE.size()-1))/fNbinsPerDecade);
	inFile.close();
}

double GZK::SecDistribution::N(int deltaE, double aK) const
{
	double kBin = log10(aK/fKmin)*fNbinsPerDecade;
	int k1 = (int)floor(kBin);
	if(kBin<0)
		return 0.;
	if(k1>=(int)fMinDeltaE.size()-1)//we are using constant sigma approximation above Kmax
	{
		return Nbin(deltaE, fMinDeltaE.size()-1);
	}
	double n1 = Nbin(deltaE, k1);
	double n2 = Nbin(deltaE, k1+1);
	double result = n1 + (n2-n1)*(kBin-(double)k1);
	ASSERT_VALID_NO(result);
	return result;
}

double GZK::SecDistribution::Nbin(int deltaE, int kBin) const
{
	const std::vector<double>& secondaries = fData(kBin);
	int i = deltaE - fMinDeltaE[kBin];
	if(i<0 || i>=(int)secondaries.size())
		return 0.;
	else
		return secondaries[i];
}

double GZK::SecDistribution::CalculateRate(const CBackgroundIntegral& aBackgr, const IFunction& aTotCs, int deltaE, double aGamma) const
{
	int kMaxBin = (int)floor(log10(2.*aBackgr.Kmax()*aGamma/fKmin)*fNbinsPerDecade);
	double sum = 0.;
	if(kMaxBin>=(int)fMinDeltaE.size())//we are using constant N approximation above Kmax
	{
		double n = Nbin(deltaE, fMinDeltaE.size()-1);
		if(n>0)
		{
			double k = fKmax*fStep;
			for(int iK=fMinDeltaE.size(); iK<=kMaxBin; iK++, k*=fStep)
				sum += n*k*k*aTotCs(k)*aBackgr.integral(0.5*k/aGamma);
		}
		kMaxBin = fMinDeltaE.size()-1;
	}

	double k = fKmin;
	for(int iK=0; iK<=kMaxBin; iK++, k*=fStep)
	{
		double n = Nbin(deltaE, iK);
		if(n>0)
			sum += n*k*k*aTotCs(k)*aBackgr.integral(0.5*k/aGamma);
	}
	sum *= (1.15129254649702/fNbinsPerDecade/aGamma/aGamma);// log(10)/fNbinsPerDecade/2/aGamma^2
	ASSERT_VALID_NO(sum);
	return sum;
}

}//namespace couplings{
