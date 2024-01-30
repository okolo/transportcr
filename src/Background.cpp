#include "Addfunc.h"
#include "Background.h"
#include <math.h>
#include "Addfunc.h"
#include "PrimackIROSpectrum.h"
#include "KneiskeIROSpectrum.h"
#include "TableFunc.h"
#include "Stecker98IROSpectrum.h"
#include "Stecker2005IROSpectrum.h"
#include "Kneiske1001IROSpectrum.h"
#include "Kneiske0309IROSpectrum.h"
#include "Franceschini08EBL.h"
#include "Sarkar1005Background.h"
#include "TimeZ.h"
#include "Units.h"
#include "Inoue12IROSpectrum.h"
#include "InjectionSpectra.h"
#include "RadioBackground.h"
// background photons spectrum

///////////////////////////////////////////////////////////////////////


/*************  neutrino background features  *****************************/

double fNeutrinoClustering=1;//neutrino clustering factor (typical 100-300). Setting to 1 switches off the clustering effect
double lNeutrinoClustering=0;//here in Mpc, further in internal units (typical 5). Setting to 0 switches off the clustering effect

/**************************************************************************/


CBackgroundIntegral::CBackgroundIntegral(CBackgroundTable& aBackground):
fTable(NULL),
fK(NULL),
fI(NULL),
fBackground(aBackground)
{
	fBackground.AddObserver(this);
}

CBackgroundIntegral::~CBackgroundIntegral()
{
	reset();
}

void CBackgroundIntegral::reset(){
	delete fTable;
	delete fK;
	delete fI;
	fTable = NULL;
	fK = NULL;
	fI = NULL;
}

//recalculate table function
void CBackgroundIntegral::OnBackgroundChange()
{
  const int nn = Ranges().nE();
  const double bmin = fBackground.Kmin();

	reset();
	CVarArray<double> x(nn);
	CVarArray<double> y(nn);
	
	fTotIntegral = 0.;

	const IFunction* F = fBackground.getBackground();
	double ac = 100.;
	double step = pow(BC().ss1,1./ac);
	double mult = log(step)/BC().sF;
	double b = fBackground.Kmax()/BC().ss1;
	for(int i=0; b>=bmin; b/=step,i++){
		fTotIntegral += F->f(b)*mult/b/b;
		//cerr << b << "\t" << F(b) << "\t" << fTotIntegral << "\n";
		x.add(b);
		y.add(fTotIntegral);
	}

	fK = new CVector(x);
	fI = new CVector(y);
	fK->invert();
	fI->invert();
	fTable = new CDefaultTableFunc(*fK,*fI);
}

double CBackgroundIntegral::calculateR(double aGamma, IFunction* aSigma, double aMinK, double aMaxK, int aAc) const
{
	if (0.5*aMinK/aGamma>=fBackground.Kmax()) return 0.;
	TBackgroundFunctional backgrI(*aSigma,*fTable,aGamma);
	double result = FuncUtils::integrate(aMinK,aMaxK,&backgrI,&gUnitFunc,aAc);
	result /= (2.*aGamma*aGamma);
	return result;
}

double CBackgroundIntegral::calculateR(double aGamma, IFunction* aSigma, double aMinK, double aMaxK, double epsrel, double epsabs, int key) const
{
	if (0.5*aMinK/aGamma>=fBackground.Kmax()) return 0.;
	TBackgroundFunctional backgrI(*aSigma,*fTable,aGamma);
	size_t maxSubIntervals = (size_t)(1./epsrel) + 1000;
	double result = funcUtils.gsl_integration_qag(&backgrI, aMinK, aMaxK, epsabs, epsrel, maxSubIntervals, key);
	result /= (2.*aGamma*aGamma);
	ASSERT_VALID_NO(result);
	return result;
}

double CBackgroundIntegral::calculateRlogsc(double aGamma, IFunction* aSigma, double aMinK, double aMaxK, double epsrel, double epsabs, int key) const
{
	if (0.5*aMinK/aGamma>=fBackground.Kmax()) return 0.;
	TBackgroundFunctionalLogsc backgrI(*aSigma,*fTable,aGamma);
	size_t maxSubIntervals = (size_t)(1./epsrel) + 1000;
	double result = funcUtils.gsl_integration_qag(&backgrI, log(aMinK), log(aMaxK), epsabs, epsrel, maxSubIntervals, key);
	result /= (2.*aGamma*aGamma);
	ASSERT_VALID_NO(result);
	return result;
}

/*
Integrate[n(b)/b^2,{b,E,Infinity}]
*/
double CBackgroundIntegral::integral(double aE) const
{
	if (aE<fBackground.Kmin()) return fTotIntegral;
	if (aE>=fBackground.Kmax()) return 0;
	return fTable->f(aE);
}

double CBackgroundIntegral::TBackgroundFunctional::f(double _x) const
{
	return _x*fBackgroundIntegral.f(0.5*_x/fGamma)*fSigma(_x);
}

double CBackgroundIntegral::TBackgroundFunctionalLogsc::f(double _lnK) const
{
	double k = exp(_lnK);
	return k*k*fBackgroundIntegral.f(0.5*k/fGamma)*fSigma(k);
}

CBackgroundTable::CBackgroundTable(IBackgroundSpectrum* aSpectrum):
fSpectrum(aSpectrum),
m_Ffit(0),
m_curZ(-1)
{
	m_Kmin = 1e300;
	m_Kmax = -1.;
	const CVector& zVals = Ranges().Z();
	for(int iZ = 0; iZ<zVals.length(); iZ++)
	{
		double minEz = aSpectrum->MinE(zVals[iZ]);
		double maxEz = aSpectrum->MaxE(zVals[iZ]);
		if(minEz<m_Kmin)
			m_Kmin = minEz;
		if(maxEz>m_Kmax)
			m_Kmax = maxEz;
	}
	//converting to internal units
	m_Kmin *= (1e-6/units.Eunit);
	m_Kmax *= (1e-6/units.Eunit);
	m_nK=Ranges().adjustRanges(m_Kmin,m_Kmax);   //.nK();
	m_backgroundEnergies.Create(m_Kmin, m_nK, BC());
	m_midBackgroundEnergies.Create(m_Kmin*BC().ss2, m_nK, BC());
	m_F.create(m_nK);
	m_Fmid.create(m_nK);
	m_F_all.create(2*m_nK);
	m_E_all.create(2*m_nK);
	for(int i=0; i<m_nK; i++)
	{
		m_E_all[2*i] = m_backgroundEnergies[i];
		m_E_all[2*i+1] = m_midBackgroundEnergies[i];
	}
}

CBackgroundTable::~CBackgroundTable()
{
	delete m_Ffit;
}

const CLinearFunc* CBackgroundTable::getBackground()
{
	if (m_Ffit) return m_Ffit;
	int maxI=nK();
	for(int i=0; i<maxI; i++)
	{
		m_F_all[2*i] = F(i);
		m_F_all[2*i+1] = Fmid(i);
	}
	m_Ffit = new CLinearFunc(m_E_all,m_F_all);
	return m_Ffit;
}


void CBackgroundTable::AddObserver(IBackgroundObserver* aObserver)
{
	m_Observers.add(aObserver);
}

void CBackgroundTable::NotifyObservers()
{
	for(int i=0; i<m_Observers.length(); i++)
		m_Observers[i]->OnBackgroundChange();
}

void CBackgroundTable::saveBackground(const char* dir) const
{
	const int mm = nK();
	char path[1024];
	char dirPath[1024];

	strcpy(dirPath,plt_local_c);

	if(dir!=NULL)
		strcat(dirPath,dir);

	if(dirPath[strlen(path)-1]==DIR_DELIMITER_CH)
		dirPath[strlen(path)-1] = '\0';

	LOG_MESSAGE2("Saving background spectrum to ", dirPath);

	char suffix[32];
	if(redshift.z()==0.)
		suffix[0] = '\0';
	else
		sprintf(suffix,"_z=%lg",redshift.z());

	sprintf(path,"%s"DIR_DELIMITER_STR"backgr%s",dirPath,suffix);
	FILE* ffspec=fopen(path,"wt");
	if(ffspec==NULL)
	{//try to create dir
		Mkdir(dirPath);
		ffspec=Fopen(path,"wt");
	}

	for(int k=0; k<mm; k++)
	{
		double b=midBackgroundE()[k];
		double Fb=Fmid(k)/BC().sF;
		if(Fb>0)
			fprintf(ffspec,"%lg  %lg\n",b/units.eV,Fb*units.cm3);
	};
	fclose(ffspec);

	sprintf(path,"%s"DIR_DELIMITER_STR"backgr.astro%s",dirPath,suffix);
	FILE* ffspec1=Fopen(path,"wt");
	fprintf(ffspec1,"# frequency (GHz)\tflux (W*m^-2*Hz^-1*sr^-1)\n");
	for(int k=0; k<mm; k++)
	{
		double b=midBackgroundE()[k];
		double Fb=Fmid(k)/BC().sF;
		if(Fb>0)
			fprintf(ffspec1,"%lg  %lg\n",b/units.Hz_photon*1e-9,Fb/4./M_PI/units.W*units.Hz_photon*units.cm*units.cm*1e4);
	};
	fclose(ffspec1);
	LOG_MESSAGE("Done");
}

void CBackgroundTable::update(double aZ/*redshift*/)
{
	if (m_curZ>=0 && aZ==m_curZ) {
		return;//no need to update
	}
	double dl = 1.+aZ;
	double xmult = units.Eunit*1e6;
	double fmult = BC().sF*units.Vunit*dl*dl*dl;
	const int mm = m_F.length();
	for(int i=0; i<mm; i++){
		m_F[i] = fSpectrum->F(m_backgroundEnergies[i]*xmult, aZ)*fmult;
		m_Fmid[i] = fSpectrum->F(m_midBackgroundEnergies[i]*xmult, aZ)*fmult;
	}
	m_curZ=aZ;

	delete m_Ffit;
	m_Ffit = 0;
	NotifyObservers();
}

CompoundBackground::CompoundBackground(double aEmin, double aEmax):
        m_Emin(aEmin),
        m_Emax(aEmax)
{
}

CompoundBackground::CompoundBackground(CompoundBackground* aCopyFrom):
m_Emin(aCopyFrom->m_Emin),
m_Emax(aCopyFrom->m_Emax)
{
    m_weights.assign(aCopyFrom->m_weights.begin(), aCopyFrom->m_weights.end());

    int nC = aCopyFrom->m_components.length();
    for (int i = 0; i < nC; i++)
        m_components.add(aCopyFrom->m_components[i]);

    m_overriden_backgr = aCopyFrom->m_overriden_backgr;
}

void CompoundBackground::addComponent(IBackgroundSpectrum* aComponent, double aWeight)
{
	ASSERT(aWeight>0);
	m_components.add(aComponent);
	m_weights.push_back(aWeight);
}

void CompoundBackground::replacePart(IBackgroundSpectrum* aComponent)
{
	m_overriden_backgr = aComponent;
}

bool CompoundBackground::init()
{
	int nC = m_components.length();
	bool result = true;
	for(int i=0; i<nC; i++)
		result = m_components[i]->init() && result;
	if((IBackgroundSpectrum*)m_overriden_backgr)
		result = m_overriden_backgr->init() && result;

    VERIFY_VALID_NO(MaxE(0) - MinE(0));

	return result;
}

double CompoundBackground::F(double E, double z)
{
    if (E < m_Emin || E>m_Emax)
        return 0;
    double result = 0.;
    if((IBackgroundSpectrum*)m_overriden_backgr && m_overriden_backgr->MaxZ()>=z && m_overriden_backgr->MinE(z)<=E &&
            m_overriden_backgr->MaxE(z)>=E)
        result = m_overriden_backgr->F(E, z);
    else {
        int nC = m_components.length();
        for (int i = 0; i < nC; i++) {
            IBackgroundSpectrum &backgr = m_components(i);
            if (backgr.MaxZ() >= z && backgr.MaxE(z) >= E && backgr.MinE(z) <= E) {
                double f = backgr.F(E, z);
                ASSERT_VALID_NO(f);
                result += f * m_weights[i];
            }
        }
    }
	ASSERT_VALID_NO(result);
	return result;
}

double CompoundBackground::MaxZ() const
{
	int nC = m_components.length();
	double result = -1;
	for(int i=0; i<nC; i++)
	{
		double maxZ = m_components[i]->MaxZ();
		if(maxZ>result)
			result = maxZ;
	}
	return result;
}

double CompoundBackground::MaxE(double aZmax) const
{
	double result = -1;
	int nC = m_components.length();
	for(int i=0; i<nC; i++)
	{
		double maxE = m_components[i]->MaxE(aZmax);
		if(maxE>result)
			result = maxE;
	}
	return m_Emax < result ? m_Emax : result;
}

double CompoundBackground::MinE(double aZmax) const
{
	double result = 1e300;
	int nC = m_components.length();
	for(int i=0; i<nC; i++)
	{
		double minE = m_components[i]->MinE(aZmax);
		if(minE<result)
			result = minE;
	}
	return m_Emin > result ? m_Emin : result;
}

ostream& operator<<(ostream& target, const CBackgroundTable& toWrite)
{
	int iMax=toWrite.nK();
	double multX=units.Eunit*1e6;// converting to eV
	double multY=1./(BC().sF*units.Vunit);// converting to n(E)*E in cm^-3
	for(int i=0; i<iMax; i++){
		if (toWrite.F(i)>0) 
			target << toWrite.backgroundE()[i]*multX << "\t" << toWrite.F(i)*multY << "\n";
		if (toWrite.Fmid(i)>0) 
			target << toWrite.midBackgroundE()[i]*multX << "\t" << toWrite.Fmid(i)*multY << "\n";
	}
	target.flush();
	return target;
}

HighRedshiftBackgrExtension::HighRedshiftBackgrExtension(IBackgroundSpectrum* aBackground, double aDeltaZconst, double aPowerLow, double aDeltaZexp, double aZmax):
		fBackground(aBackground),
		fPowerLow(aPowerLow),
		fDeltaZexp(aDeltaZexp),
		fZmax(aZmax),
		fInnerZmax(0)
{
	fZconst=aDeltaZconst;
}

bool HighRedshiftBackgrExtension::init()
{
	if(fInnerZmax>0)
		return true;//already initialized
	if(!fBackground->init())
		return false;
	fInnerZmax = fBackground->MaxZ();
	fZconst+=fInnerZmax;
	return true;
}

double HighRedshiftBackgrExtension::MaxZ() const
{
	return fZmax;
}

double HighRedshiftBackgrExtension::MaxE(double aZmax) const
{
	return fBackground->MaxE(MIN(aZmax,fZmax));
}

double HighRedshiftBackgrExtension::MinE(double aZmax) const
{
	return fBackground->MinE(MIN(aZmax,fZmax));
}

double HighRedshiftBackgrExtension::F(double aE, double aZ)
{
	if(aZ<fInnerZmax)
		return fBackground->F(aE, aZ);
	if(aZ>fZmax)
		return 0.;
	if(aZ<fZconst)
		return fBackground->F(aE, 0.9999*fInnerZmax);

	double powerLowFactor = pow((1.+aZ)/(1.+fZconst),fPowerLow);
	double expFactor = (fDeltaZexp>0)?exp((fZconst-aZ)/fDeltaZexp):1;
	return powerLowFactor*expFactor*fBackground->F(aE, 0.9999*fInnerZmax);
}

CModifiedBlackbodySpectrum::CModifiedBlackbodySpectrum(double aTemperature/*[T]=K*/, double aPower, double aDensity, bool aZdependence, double aZmax, double aWidthDecade):
iZmax(aZmax),
iPower(aPower),
iMultiplier(1.),
iT(aTemperature/1.16e10/units.Eunit),//internal units
iZdependence(aZdependence),
iRightWidthMult(pow(10., aWidthDecade/3.)),
iLeftWidthMult(pow(10., 2.*aWidthDecade/3.)),
iDensity(aDensity)
{
}

bool CModifiedBlackbodySpectrum::init()
{
	if(iDensity>0)
	{//normalize at z=0
		iMultiplier = 1.;
		double n = BackgroundUtils::CalcIntegralDensity(*this,0);
		iMultiplier = iDensity/n;
	}
	return true;
}

double CModifiedBlackbodySpectrum::MaxE(double aZmax) const
{
	double T = iZdependence ? iT*(1.+aZmax) : iT;
	double maxK = iRightWidthMult<700 ? iRightWidthMult*T : 700*T; //limited by double accuracy (see method F(E,z))
	return maxK*units.Eunit*1e6;//converting to eV
}

double CModifiedBlackbodySpectrum::MinE(double aZmax) const
{
	return iT/iLeftWidthMult*units.Eunit*1e6;
}

double CModifiedBlackbodySpectrum::F(double E/*eV*/, double z)
{
	E = E/1e6/units.Eunit;
	double dl = (1.+z);
	double gamma=iZdependence?E/iT/dl:E/iT;
	double result = 0.;
	if(gamma<700)//limited by double accuracy
	{
		double mult = iPower==0. ? iMultiplier : iMultiplier*pow(gamma, iPower);
		double ex_1 = (gamma > 1e-2) ? Exp(gamma)-1.0 : gamma*(1.+ 0.5*gamma*(1. + 0.333333333333333*gamma*(1. + 0.25*gamma)));
		result = mult*E*E*E/9.869604401089/ex_1;//that is E^3/Pi^2/(...) - two internal degrees of freedom assumed
	}
	result /= (units.Vunit*dl*dl*dl);//dividing by dv to make it in comoving frame
	return result;
}

double BackgroundUtils::CalcIntegralDensity(IBackgroundSpectrum& aBackground, double aZ, double aRelError)
{
	class Kern : public IFunction
	{
		double Z;
		IBackgroundSpectrum& Background;
	public:
		Kern(IBackgroundSpectrum& aBackground, double aZ):Z(aZ),Background(aBackground){}
		virtual double f(double aLogE) const
		{
			double E = exp(aLogE);
			double result = Background.F(E, Z);
			ASSERT_VALID_NO(result);
			return result;
		}
	};
	Kern k(aBackground, aZ);
	double comovDensity = funcUtils.gsl_integration_qag(&k, log(aBackground.MinE(aZ)), log(aBackground.MaxE(aZ)), 0, aRelError, 100+floor(1./aRelError));
	double dl = 1.+aZ;
	return comovDensity*dl*dl*dl;
}

double BackgroundUtils::CalcIntegralPowerDensity(IBackgroundSpectrum& aBackground, double aZ, double aRelError)
{
	class Kern : public IFunction
	{
		double Z;
		IBackgroundSpectrum& Background;
	public:
		Kern(IBackgroundSpectrum& aBackground, double aZ):Z(aZ),Background(aBackground){}
		virtual double f(double aLogE) const
		{
			double E = exp(aLogE);
			return Background.F(E, Z)*E;
		}
	};
	Kern k(aBackground, aZ);
	double comovDensity = funcUtils.gsl_integration_qag(&k, log(aBackground.MinE(aZ)), log(aBackground.MaxE(aZ)), 0, aRelError, 100+floor(1./aRelError));
	double dl = 1.+aZ;
	return comovDensity*dl*dl*dl;
}
