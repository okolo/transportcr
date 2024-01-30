/*
 * AA.cpp
 *
 *  A-A interaction (experimental, work in progress)
 */

#include <cmath>
#include "AA.h"
#include "Medium.h"
#include "Units.h"
#include "Nucleus.h"
#include "TableReader.h"

#ifndef NO_FORTRAN
extern "C"
{
double spec_gam_(double& E_p, double& E_g, int& reac/*=1*/);
double spec_el_qgs_(double& E_p,double& E_el,int& reac/*=1*/);
double spec_nu_qgs_(double& E_p,double& E_nu,int& reac/*=1*/,int& fl);
void pp_init_();
}
#else
double spec_gam_(double& E_p, double& E_g, int& reac/*=1*/)
{
	ThrowError("Fatal error in spec_gam_: Fortran code is disabled");
	return 0.;
}

double spec_el_qgs_(double& E_p,double& E_el,int& reac/*=1*/)
{
	ThrowError("Fatal error in spec_el_qgs_: Fortran code is disabled");
	return 0.;
}

double spec_nu_qgs_(double& E_p,double& E_nu,int& reac/*=1*/,int& fl)
{
	ThrowError("Fatal error in spec_nu_qgs_: Fortran code is disabled");
	return 0.;
}

void pp_init_()
{
	ThrowError("Fatal error in pp_init_: Fortran code is disabled");
}
#endif

namespace couplings {

#define PP_CS_DATA_DIR DATA_DIR "pp/cross_sec/"

void NucleonProton::InitWithTableRate()
{
	if(fTableFile.length()==0)
		ThrowError("p-p rate table is not set");
	const int nn = Ranges().nE();
	CTableReader tr(fTableFile,2);
	CVector E(tr.getColumn(0));
	CVector tau(tr.getColumn(1));
	E *= (1e-6/units.Eunit);//converting eV to internal units
	tau *= (units.Lunit/units.Mpc_cm);//converting tau for 1 Mpc to rate in internal units
	CLogScaleLinearFunc tauF(E, tau);
	int iE = Ranges().nMinInteractionN();
	if(E[0]>Ranges().midE()[iE] && E[0]>fEth)
	{//extrapolating tau to lower energies assuming power law dependence of grammage on energy
		double X0 = tau[0]/Sigma_inel(E[0]);//calculating grammage divided by m_p
		double X1 = tau[1]/Sigma_inel(E[1]);
		double alpha = log(X0/X1)/log(E[0]/E[1]);
		double C = X0/pow(E[0],alpha);
		for(;Ranges().midE()[iE]<=E[0]; iE++)
		{
			double curE = Ranges().midE()[iE];
			double X = C * pow(curE, alpha);
			double curSigma = Sigma_inel(curE);
			fRates[iE] = curSigma*X;
		}
	}
	////extrapolating above table Emax assuming power law dependence of grammage on energy
	double E1 = E[E.length()-1];
	double E0 = E[E.length()-2];
	double X1 = tau[E.length()-1]/Sigma_inel(E1);
	double X0 = tau[E.length()-2]/Sigma_inel(E0);
	double alpha = log(X0/X1)/log(E0/E1);
	double C = X0/pow(E0,alpha);
	for(; iE<nn; iE++)
	{
		double curE = Ranges().midE()[iE];
		if(curE>E1)
		{//extrapolating assuming power law dependence of grammage on energy
			double X = C * pow(curE, alpha);
			double curSigma = Sigma_inel(curE);
			fRates[iE] = curSigma*X;
		}
		else
			fRates[iE] = tauF(curE);
	}
}

void NucleonProton::InitWithTableGrammage()
{
	double m = ParticleData::getParticleMass(EProton);
	if(fTableFile.length()==0)
		ThrowError("p-p grammage table is not set");
	const int nn = Ranges().nE();
	CTableReader tr(fTableFile,2);
	CVector E(tr.getColumn(0));
	CVector N(tr.getColumn(1));
	E *= (1e-6/units.Eunit);//converting eV to internal units
	N *= (units.Vunit*units.gram/m/units.Mpc_cm);//converting 1Mpc grammage to concentration in internal units
	//for 1 Mpc to rate in internal units
	CLogScaleLinearFunc conc(E, N);
	int iE = Ranges().nMinInteractionN();
	if(E[0]>Ranges().midE()[iE] && E[0]>fEth)
	{//extrapolating grammage to lower energies assuming gramage ~ E^{-1/3}
		double alpha = -1./3.;
		double C = N[0]/pow(E[0],alpha);
		for(;Ranges().midE()[iE]<=E[0]; iE++)
		{
			double curE = Ranges().midE()[iE];
			double X = C * pow(curE, alpha);
			fRates[iE] = X*Sigma_inel(curE);
		}
	}
	////extrapolating rates above table Emax assuming constant grammage
	double maxTableE = E[E.length()-1];
	double limX = N[E.length()-1];
	for(; iE<nn; iE++)
	{
		double curE = Ranges().midE()[iE];
		if(curE>maxTableE)
		{//extrapolating assuming constant grammage
			fRates[iE] = limX*Sigma_inel(curE);
		}
		else{
			fRates[iE] = conc(curE)*Sigma_inel(curE);
		}
	}
}

void NucleonProton::SetBackgrounds(const Medium& aPropagCoef)
{
	const int nn = Ranges().nE();
	if(!fRates.length())
		fRates.create(nn);
	if(fPconc==aPropagCoef.protonsConc)
		return;
	if(fPconc>0)
	{
		fRates *= (aPropagCoef.protonsConc/fPconc);
		fPconc = aPropagCoef.protonsConc;
		return;
	}

	fPconc = aPropagCoef.protonsConc;
	if(fMode==RMdefault)
	{///standard mode without table rate
		for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
		{
			double E = Ranges().midE()[iPrim];
			fRates[iPrim] = Sigma_inel(E)*fPconc;
		}
	}
	else if(fMode==RMtau)
	{///rate table is given
		InitWithTableRate();
	}
	else if(fMode==RMgrammage)
	{///rate table is given
		InitWithTableGrammage();
	}
}

NucleonProton::NucleonProton(RateMode mode, PPCrossSectionMode aScMode, const char* aTable):
		fPconc(0),
		fEth(1.22e3/units.Eunit),
		fTableFile(aTable==0?"":aTable),
		fMode(mode),
		fScMode(aScMode),
		fTableSigmaLeftMult(0),
		fTableSigmaRightMult(0)
{
	TParticle primaries[] = {ENeutron, EProton};
	TParticle secondaries[] = {EPhoton, ENeutrinoM, ENeutrinoAM, ENeutrinoE, ENeutrinoAE, EElectron, EPositron, EEndParticle};
	for(int p=0; p<2; p++)
	{
		TParticle primary = primaries[p];
		if(aTable && primary == ENeutron )
			continue;//use table p-p for protons only
		AddChannel(new Channel_N_N(this, primary));
		if(aScMode == PPAnalKelner)
		{
			AddChannel(new Channel_N_gamma(this, primary));
			AddChannel(new Channel_N_mnu(this, primary, ENeutrinoM));
			AddChannel(new Channel_N_mnu(this, primary, ENeutrinoAM));
			AddChannel(new Channel_N_e_ne(this, primary, ENeutrinoE));
			AddChannel(new Channel_N_e_ne(this, primary, ENeutrinoAE));
			AddChannel(new Channel_N_e_ne(this, primary, EElectron));
			AddChannel(new Channel_N_e_ne(this, primary, EPositron));
		}
		else
		{
			fTableSigma = new CLogScaleLinearFunc(new CTableReader(PP_CS_DATA_DIR "sigma_in", 2));
			//use analytical cross sec dependence on E outside table range but make sure sigma is continuous
			if(fTableSigma->firstX()>fEth/units.GeV)
				fTableSigmaLeftMult = (fTableSigma->firstF()*units.mbarn)/Sigma_inel_anal(fTableSigma->firstX()*units.GeV);
			fTableSigmaRightMult = (fTableSigma->lastF()*units.mbarn)/Sigma_inel_anal(fTableSigma->lastX()*units.GeV);

			//fTableSigmaRightMult=1; fTableSigmaLeftMult=1;//test
			for(int sec=0; secondaries[sec]<EEndParticle; sec++)
			{
				if(aScMode == PPTableKachelriess)
				{
					ThrowError("not implemented:\nTODO: prepare p-p tables by converting original M.Kachelriess tables to matrix format");
					AddChannel(new Channel_N_xTable(this, primary, secondaries[sec]));
				}
				else if(aScMode == PPCodeKachelriess)
					AddChannel(new Channel_N_xKachelriess(this, primary, secondaries[sec]));
				else
					ThrowError("pp: unsupported cross section mode");
			}
		}
	}
}

bool NucleonProton::Channel_N_xKachelriess::fInitialized=false;

NucleonProton::Channel_N_xKachelriess::Channel_N_xKachelriess(NucleonProton* aCoupling,
		TParticle aPrimary, TParticle aSecondary):
		CouplingChannelT<NucleonProton>(aCoupling, aPrimary, aSecondary)
{
	if(!fInitialized)
	{
		pp_init_();
		fInitialized = true;
	}
}

void NucleonProton::Channel_N_xKachelriess::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	double tableEmin = fTableEmin_eV[Secondary()]*units.eV;
	double tableEmax = fTableEmax_eV[Secondary()]*units.eV;
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{//TODO: optimize (build table to avoid coefficient recalculation)
		double rate = fCoupling.fRates[iPrim];
		if(rate<=0)
			continue;
		ASSERT_VALID_NO(rate);
		rate *= (BC().sF);
		double E = Ranges().midE()[iPrim];
		double sigma_tot = 0;
		if(E<tableEmin)
		{//extrapolation will be assuming dSigma/dE_sec * E_sec / sigma_tot independent on E_p
			sigma_tot = fCoupling.Sigma_inel(tableEmin);
		}
		else if(E>tableEmax)
		{//extrapolation will be assuming dSigma/dE_sec * E_sec / sigma_tot independent on E_p
			sigma_tot = fCoupling.Sigma_inel(tableEmax);
		}
		else
		{
			sigma_tot = fCoupling.Sigma_inel(E);
		}
		if(sigma_tot<=0)
			ThrowError("pp: sigma==0 while rate > 0");
		rate /= sigma_tot;

		for(int iSec = 0; iSec<=iPrim; iSec++)
		{
			double Esec = Ranges().midE()[iSec];
			double N = 0.;
			if(E<tableEmin)
			{
				N = sigma(tableEmin, Esec*tableEmin/E);
			}
			else if(E>tableEmax)
			{
				N = sigma(tableEmax, Esec*tableEmax/E);
			}
			else
			{
				N = sigma(E, Esec);
			}
			ASSERT_VALID_NO(N);
			if(N>0)
				aCoef.Add(iSec, iPrim, N*rate);
		}
	}
}

//range of proton energies in differential cross section tables
const double NucleonProton::fTableEmax_eV[] = {2.5e18,2.5e18,2.5e18,2.5e18,2.5e18,-1,2.5e18,2.5e18};
const double NucleonProton::fTableEmin_eV[] = {15848000049.591064,15848000049.591064,15848000049.591064,1.58480005E+10,1.58480005E+10,-1,1.58480005E+10,1.58480005E+10};

double NucleonProton::Channel_N_xKachelriess::sigma(double Eprim, double Esec) const
{
	Eprim /= units.eV;
	Esec /= units.eV;
	int reac = 1;
	int fl = 1;
	switch(Secondary())
	{
	case EPhoton:
	{
		RERURN_VALID(spec_gam_(Eprim, Esec, reac)*units.mbarn);
	}
	case EElectron:
	case EPositron:
	{
		RERURN_VALID(0.5*spec_el_qgs_(Eprim, Esec, reac)*units.mbarn);//table was built for sigma(electrons) + sigma(positrons)
	}
	case ENeutrinoE:
		fl = 1;
		break;
	case ENeutrinoAE:
		fl = 2;
		break;
	case ENeutrinoM:
		fl = 3;
		break;
	case ENeutrinoAM:
		fl = 4;
		break;
	default:
		ASSERT(false);
		return 0.;
	}
	RERURN_VALID(spec_nu_qgs_(Eprim, Esec, reac, fl)*units.mbarn);
}

NucleonProton::Channel_N_xTable::Channel_N_xTable(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary):
		CouplingChannelT<NucleonProton>(aCoupling, aPrimary, aSecondary),
		fSigma(ToString(PP_CS_DATA_DIR)+ParticleData::getParticleFileName(aSecondary))

{
}

void NucleonProton::Channel_N_xTable::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
	double diffSigmaUnits = units.mbarn;
	double logStep = 1./BC().s;
	double logEmin_eV = log10(Ranges().midE()[0]*units.Eunit*1e6);
	double logE_eV = logEmin_eV+logStep*Ranges().nMinInteractionN();
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++, logE_eV+=logStep)
	{//TODO: optimize (build table to avoid coefficient recalculation)
		double rate = fCoupling.fRates[iPrim];
		if(rate<=0)
			continue;

		rate *= (diffSigmaUnits*BC().sF);
		double tableSigma = fCoupling.Sigma_inel_tab(Ranges().midE()[iPrim]);
		double logEsec_eV = logEmin_eV;

		for(int iSec = 0; iSec<=iPrim; iSec++, logEsec_eV+=logStep)
		{
			double N = fSigma(logE_eV,logEsec_eV)/tableSigma;
			if(N>0)
				aCoef.Add(iSec, iPrim, N*rate);
		}
	}
}


double NucleonProton::Sigma_inel(double aEp) const
{
	if(aEp<fEth)
		return 0.;
	if(fScMode == PPAnalKelner)
		return Sigma_inel_anal(aEp);
	double result = Sigma_inel_tab(aEp);
	if(result==0 && aEp>fEth)
		return Sigma_inel_anal(aEp);//TODO: invent something smarter
	else
		return result;
}

double NucleonProton::Sigma_inel_anal(double aEp) const
{//formula (79) of astro-ph/0606058v1
	if(aEp<fEth)
		return 0.;
	double L = log(aEp/units.TeV);
	double r = fEth/aEp;
	r*=r;//(aEp/fEth)^2
	double result = 34.3+L*(1.88+0.25*L)*sqrt(1.-r*r);
	result *= units.mbarn;
	ASSERT_VALID_NO(result);
	return result;
}

double NucleonProton::Sigma_inel_tab(double aEp) const
{
	if(aEp<fEth)
		return 0.;
	double Ep_GeV = aEp/units.GeV;//table E is in GeV
	double result = 0;
	if(Ep_GeV<fTableSigma->firstX())
	{//use analytical cross sec dependence on E but make sure sigma is continuous
		result = fTableSigmaLeftMult*Sigma_inel_anal(aEp);
	}
	else if(Ep_GeV>fTableSigma->lastX())
	{//use analytical cross sec dependence on E but make sure sigma is continuous
		result = fTableSigmaRightMult*Sigma_inel_anal(aEp);
	}
	else
	{
		result = (*fTableSigma)(Ep_GeV)*units.mbarn;//table sigma is in mb
	}
	ASSERT_VALID_NO(result);
	return result;
}

void NucleonProton::Channel_N_N::Coef(CMatrixAddOnlyView& aCoef) const
{
	double n_gas = fCoupling.fPconc;
	if(n_gas<=0)
		return;
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		double rate = fCoupling.fRates[iPrim];
		if(rate>0)
			aCoef.Add(iPrim, iPrim, -rate);
	}
}

void NucleonProton::Channel_N_x::Coef(CMatrixAddOnlyView& aCoef) const
{
	double n_gas = fCoupling.fPconc;
	if(n_gas<=0)
		return;
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{//TODO: optimize (build table to avoid coefficient recalculation)
		double E = Ranges().midE()[iPrim];
		double rate = fCoupling.fRates[iPrim];
		if(rate<=0)
			continue;
		rate *= (BC().sF);
		double L=log(E*units.Eunit*1e-6);//log(E_p/TeV)
		for(int iSec = 0; iSec<=iPrim; iSec++)
		{
			double x = Ranges().midE()[iSec]/E;
			double f = F(x,L);
			if(f>0)
				aCoef.Add(iSec, iPrim, x*rate*f);
		}
	}
}

double NucleonProton::F_gamma(double x, double L)
{
	if(x>=1.)
		return 0.;
	double B=1.30+0.14*L+0.011*L*L;
	double beta = 1./(1.79+0.11*L+0.008*L*L);
	double k = 1./(0.801+0.049*L+0.014*L*L);
	double lnx = log(x);
	double xb = pow(x,beta);

	double A = (1.-xb)/(1.+k*xb*(1.-xb));
	A*=A;

	double result = B*lnx/x*A*A*(1./lnx-4.*beta*xb*(1./(1.-xb)+k*(1.-2.*xb)/(1.+k*xb*(1.-xb))));
	ASSERT_VALID_NO(result);
	return result;
}

double NucleonProton::F_mnu1(double x, double L)
{
	double y = x/0.427;
	if(y>=1)
		return 0.;
	double lny = log(y);
	double B = 1.75 + 0.204*L + 0.01*L*L;
	double beta = 1./(1.67+0.111*L+0.0038*L*L);
	double k = 1.07 - 0.086*L + 0.002*L*L;
	double yb = pow(y,beta);

	double A = (1.-yb)/(1.+k*yb*(1.-yb));
	A*=A;

	double result = B*lny/y*A*A*(1./lny-4.*beta*yb*(1./(1.-yb)+k*(1.-2.*yb)/(1.+k*yb*(1.-yb))));
	ASSERT_VALID_NO(result);
	return result;
}

double NucleonProton::F_e(double x, double L)
{
	if(x>=1.)
		return 0.;
	double B = 1./(69.5+2.65*L+0.3*L*L);
	double b = 0.201+0.062*L+0.00042*L*L;
	if(b<0)
		return 0.;
	double beta = pow(b,-0.25);
	double L23 = 2.3+L;
	double k = (0.279+0.141*L+0.0172*L*L)/(0.3+L23*L23);
	double lnx = log(x);
	double lnx2 = lnx*lnx;
	double xb = pow(x,beta);
	double A = 1.+k*lnx2;
	A *= (A*A);

	double result = -B*A/(x*(1.+0.3/xb))*lnx2*lnx2*lnx;
	ASSERT_VALID_NO(result);
	return result;
}

#define AP_CS_DATA_DIR DATA_DIR "AA/"

	NucleiProton::NucleiProton(const char *aGrammageTable) : fRateMult(1.), fLastConc(1.){
		if(aGrammageTable && strlen(aGrammageTable)==0)
			aGrammageTable = 0;
		const double pi0mass = 134.9766*units.MeV;
		const double protonM = ParticleData::getParticleMass(EProton);
		CParticleList plist;
		int maxA = EEndNuclei-EStartNuclei+1;

		CAutoDeletePtrArray<CVector> data;
		CAutoDeletePtrArray<CLogScaleLinearFunc> sigmas;
		CVarArray<int> dataA;
		//reading data
		for(int A=0; A <= maxA ; A++) {//here A=0 correspond to neutron
			TParticle nucl = (TParticle) (ENeutron + A);
			double scaleFactorE = ParticleData::GetEnergyScaleFactor(nucl);
			std::string csFile = AP_CS_DATA_DIR + ToString(A);
			if (fileExists(csFile)) {
				CTableReader tr(csFile, 4);
				CVector *colE = new CVector(tr.getColumn(2));
				CVector *colSigma(new CVector(tr.getColumn(3)));
				data.add(colE);
				data.add(colSigma);
				(*colE) *= (units.GeV);//converting to energy in internal units
				(*colSigma) *= units.mbarn;
				CLogScaleLinearFunc* sigma = new CLogScaleLinearFunc(*colE, *colSigma);
				sigma->SetExtension(ExtConst, false, true);
				sigma->SetExtension(ExtLinear, true, false);//will extrapolate linearily in logscale till threshold energy
				sigmas.add(sigma);
				dataA.add(A);
			}
			else{
				if(A<=1 || A==maxA)//we use interpolation only for 1<A<Amax; tables for p and n and Fe must exist
					ThrowError(ToString("Missing table for p-") + ParticleData::getParticleFileName(nucl) + " interaction");
			}
		}

		//creating rate tables
		int nn = Ranges().nE();
		int nn0 = Ranges().nMinInteractionN();
		const CVector& p = Ranges().midE();//momentum
		fRatesPerUnitDensity.create(maxA+1, nn);
		//Interpolate cross section for intermediate A and E using linear interpolation in A^(2/3) and log interpolation in momentum
		int nextTableIndex = 0;
		for(int A=0; A <= maxA ; A++) {//here A=0 correspond to neutron
			TParticle nucl = (TParticle)(ENeutron + A);
			double nuclM = ParticleData::getParticleMass(nucl);
			double protonM = ParticleData::getParticleMass(EProton);
			double thresholdE = pi0mass*((0.5*pi0mass+nuclM)/protonM+1.) + nuclM;
			//double thresholdP = sqrt(thresholdE*thresholdE-nuclM*nuclM);
			double scaleFactorP = ParticleData::GetEnergyScaleFactor(nucl);
			CVector &rates = fRatesPerUnitDensity[A];
			std::string csFile = AP_CS_DATA_DIR + ToString(A);
			if (dataA[nextTableIndex]==A) {//cross section data exists
				CLogScaleLinearFunc& sigma = sigmas(nextTableIndex);
				for (int iP = nn0; iP < nn; iP++) {
					double curP = p[iP]* scaleFactorP;
					if(curP >thresholdE)
						rates[iP] = sigma(curP);//TODO: ask Michael about meaning of E column in data file
					// is it momentum or E-m (energy - mass) or energy per nucleon (it is for sure not the total energy since
					// for heavy nuclei the data also starts from 15 GeV which is less than their mass
				}
				nextTableIndex++;
			}
			else{//using linear interpolation in A^(2/3)
				double interpLaw = 2./3.;
				double x1 = pow(dataA[nextTableIndex-1],interpLaw);
				double x2 = pow(dataA[nextTableIndex],interpLaw);
				double x = pow(A,interpLaw);
				const CLogScaleLinearFunc& sigma1 = sigmas(nextTableIndex-1);
				const CLogScaleLinearFunc& sigma2 = sigmas(nextTableIndex);
				double frac1 = (x2-x)/(x2-x1);
				double frac2 = (x-x1)/(x2-x1);
				for (int iP = nn0; iP < nn; iP++) {
					double curP = p[iP]* scaleFactorP;
					if(curP >thresholdE)
						rates[iP] = frac1*sigma1(curP) + frac2*sigma2(curP);
				}
			}
		}

		//adjusting rates if grammage is given
		if(aGrammageTable) {
			std::string fTableFile = aGrammageTable;
			const int nn = Ranges().nE();
			CTableReader tr(fTableFile, 2);
			CVector colE(tr.getColumn(0));
			CVector colN(tr.getColumn(1));

			colE *= units.eV;//converting eV to internal units
			colN *= (units.Vunit * units.gram / protonM /
				  units.Mpc_cm);//converting 1Mpc grammage to concentration in internal units
			//for 1 Mpc to rate in internal units
			CLogScaleLinearFunc conc(colE, colN);
			conc.SetExtension(ExtConst, false, true);
			conc.SetExtension(ExtLinear, true, false);

			for(int A=1; A <= maxA ; A++){//we assume that grammage is function of E/Z
				TParticle nucl = (TParticle)(ENeutron + A);
				double scaleFactorE = ParticleData::GetEnergyScaleFactor(nucl);
				double Z = CNucleus::getZ(nucl);
				CVector &rates = fRatesPerUnitDensity[A];
				for (int iE = nn0; iE < nn; iE++) {
					rates[iE] *= conc(p[iE]*scaleFactorE/Z);
				}
			}
			(fRatesPerUnitDensity[0]) *= conc.lastF();//using high energy grammage limit for neutruns (Z=0)
			fLastConc = -1.;//don't multiply by Background()->protonsConc;
		}

		for(int particle=ENeutron; particle<EEndNuclei; particle++){
			AddChannel(new Channel_A_A(this,(TParticle)particle));
		}
	}

	void NucleiProton::Channel_A_A::Coef(CMatrixAddOnlyView& aCoef) const{
		if(fCoupling.fRateMult<=0)
			return;

		const int nn = Ranges().nE();
		CVector& rates = fCoupling.fRatesPerUnitDensity[Primary()-ENeutron];
		for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
		{
			aCoef.Add(iPrim, iPrim, -rates[iPrim]*fCoupling.fRateMult);
		}
	}

	void NucleiProton::SetBackgrounds(const Medium& aPropagCoef){
		double n_p = aPropagCoef.protonsConc;
		ASSERT(n_p>0);///avoid division by 0 in future
		if(fLastConc<0)//support for grammage
			fRateMult = 1.;
		else
			fRateMult = fRateMult*(n_p/fLastConc);
		fLastConc = n_p;
	}

}
