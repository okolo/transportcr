/*
 * Galaxy.cpp
 *
 *  Created on: Aug 4, 2014
 *      Author: mac
 */

#include <math.h>
#include "Galaxy.h"
#include "Addfunc.h"
#include "Units.h"
#include "TimeZ.h"
#include "Jacobian.h"
#include "TableReader.h"

//#include <omp.h>

DiskCenterCS::DiskCenterCS():
sinBeta(Z_sun/R_sun),
cosBeta(sqrt(1.-sinBeta*sinBeta))
{
}

void DiskCenterCS::GalacticToDiskCenterCS(double l, double b, double r, double &x, double &y, double &z) const
{
	//convert to radians
	l *= (M_PI/180.);
	b *= (M_PI/180.);
	//galactic rectangular coordinates
	double zGal = r*sin(b);
	double xGal = r*cos(b);
	double yGal = xGal;
	xGal *= cos(l);
	yGal *= sin(l);
	//parallel shift of CO center to Milky Way center
	xGal -= R_sun;
	//CO rotation
	x = xGal*cosBeta-zGal*sinBeta;
	y = yGal;
	z = xGal*sinBeta+zGal*cosBeta;
}

double DiskCenterCS::GetDistanceFromSun(double l, double b, double aDistanceFromDiskCenterKpc) const
{//solving quadratic equation
	//if aDistanceFromDiskCenterKpc<R_sun two solutions exist
	//the largest distance is returned by this function

	double prodCos = cos(l*M_PI/180.)*cos(b*M_PI/180.);
	// Discriminant/4
	double D = aDistanceFromDiskCenterKpc*aDistanceFromDiskCenterKpc - R_sun*R_sun*(1.-prodCos*prodCos);
	if(D<0)
		ThrowError("DiskCenterCS.GetDistanceFromSun : invalid arguments");//check geometry
	double result = R_sun*prodCos + sqrt(D);
	ASSERT_VALID_NO(result);
	return result;
}

//Calculate distance from galactic center using galactic coordinates l,b and distance from sun rKpc in kpc
double DiskCenterCS::GetDistanceFromCenter(double l, double b, double rKpc) const
{
	//convert to radians
	l *= (M_PI/180.);
	b *= (M_PI/180.);
	//galactic rectangular coordinates
	double zGal = rKpc*sin(b);
	double xGal = rKpc*cos(b);
	double yGal = xGal;
	xGal *= cos(l);
	yGal *= sin(l);
	//parallel shift of CO center to Milky Way center
	xGal -= R_sun;
	return sqrt(xGal*xGal+yGal*yGal+zGal*zGal);
}

const double DiskCenterCS::Z_sun = 2.5e-3;//in kpc
const double DiskCenterCS::R_sun = 8.5;//in kpc

enum GalDistribution
{
	Stellar = 0,
	NFW,
	HotGas,
	HotGas2,
    GDConst,
	GDBurkert,

	EndGalDistribution
};

GalacticDensityDistribution::~GalacticDensityDistribution()
{

}

class ConstGalDistribution : public GalacticDensityDistribution, DiskCenterCS
{
public:
    ConstGalDistribution(double aRmaxKpc = 100.):
            fRmaxKpc(aRmaxKpc)
    {
    }

    double Density(double l, double b, double rKpc) const
    {
        double r = GetDistanceFromCenter(l, b, rKpc);
        return r<fRmaxKpc?1.:0.;
    }

    double MaxDistanceKpc(double l, double b) const
    {
        return GetDistanceFromSun(l, b, fRmaxKpc);
    }

    bool IsIsotropic() const{
        return true;
    }

private:
    double fRmaxKpc;
};

class StellarDistribution : public GalacticDensityDistribution, DiskCenterCS
{
	inline double n_h(double rho,double z) const
	{
		const double N_h = 2.77;
		const double f_h = 0.0051;
		const double q_h = 0.64;
		return f_h*pow(R_sun/sqrt(rho*rho+z*z/q_h/q_h),N_h);
	}
	inline double n_1(double rho, double z, double L, double H) const
	{
		return exp((R_sun-rho)/L - (fabs(z)+Z_sun)/H);
	}
	inline double n_d(double rho, double z) const
	{
		const double L1 = 2.6;//kpc
		const double L2 = 3.6;//kpc
		const double H1 = 0.3;//kpc
		const double H2 = 0.9;//kpc
		const double f_f = 0.12;
		return n_1(rho, z, L1, H1) + f_f*n_1(rho, z, L2, H2);
	}
	double Density(double l, double b, double rKpc) const
	{
		double x,y,z;
		GalacticToDiskCenterCS(l, b, rKpc, x, y, z);
		double rho = sqrt(x*x+y*y);
		return n_d(rho, z) + n_h(rho,z);
	}

	double MaxDistanceKpc(double l, double b) const
	{
		const double MaxR = 15.;//kpc
		return GetDistanceFromSun(l, b, MaxR);
	}

	bool IsIsotropic() const{
		return true;
	}
};

class NFWDistribution : public GalacticDensityDistribution, public DiskCenterCS
{
public:
	NFWDistribution(double aSingularitySmoothSizeKpc):
		C_v(12.),
		R_v(260.),
		R_s(R_v/C_v),
		SingularitySmoothSize(aSingularitySmoothSizeKpc/R_s)
	{}
	double MaxDistanceKpc(double l, double b) const
	{
		return GetDistanceFromSun(l, b, R_v);
	}

	bool IsIsotropic() const{
		return true;
	}

	double Density(double l, double b, double rKpc) const
	{
		double x = GetDistanceFromCenter(l, b, rKpc)/R_s;
		if(x<SingularitySmoothSize)
			x = SingularitySmoothSize;
		return 1./x/(1.+x)/(1.+x);
	}
protected:
	const double C_v;
	const double R_v;
	const double R_s;
	const double SingularitySmoothSize;
};

class BurkertDistribution : public GalacticDensityDistribution, public DiskCenterCS
{
public:
	BurkertDistribution():
			R_s(30.28),//as if in arxiv/1012.4515v4 Fig.1
			R_v(260)
	{

	}
	double MaxDistanceKpc(double l, double b) const
	{
		return GetDistanceFromSun(l, b, R_v);
	}

	bool IsIsotropic() const{
		return true;
	}

	double Density(double l, double b, double rKpc) const
	{
		double x = GetDistanceFromCenter(l, b, rKpc)/R_s;

		return 1./(1.+x)/(1.+x*x);
	}
protected:
	const double R_s;
	const double R_v;
};

//hot gas distribution multiplied by cosmic rays distribution (1/r assumed for CR)
class HotGasDistributionTimesCR : public NFWDistribution, protected Parameters
{
public:
	HotGasDistributionTimesCR(double aSingularitySmoothSizeKpc):
		NFWDistribution(aSingularitySmoothSizeKpc),
		C_c(1.),
		fA(1.)
	{
		Reader()->readDoublePar("HG_a", fA);
		Reader()->readDoublePar("HG_Cc", C_c);
	}

	double MaxDistanceKpc(double l, double b) const
	{
		return GetDistanceFromSun(l, b, R_v);
	}

	bool IsIsotropic() const{
		return true;
	}

	double Density(double l, double b, double rKpc) const
	{
		double x = GetDistanceFromCenter(l, b, rKpc)/R_s;
		if(x<SingularitySmoothSize)
			x = SingularitySmoothSize;
		double result = 0.;
		if(x<C_c)//i.e. if R<R_c (where R_c = R_s*C_c)
			result = n1(x);
		else
			result = C_c*C_c/(x*x);//make density continuous function provided that n1(C_c)===1.
		result /= pow(x,fA); //multiplying by cosmic rays density
		ASSERT_VALID_NO(result);
		return result;
	}
private:
	double n1(double x) const
	{
		return pow(1.+3.7/x*log(1.+x)-3.7/C_c*log(1.+C_c), 1.5);
	}
	double C_c;
	double fA;
};

//hot gas distribution multiplied by cosmic rays distribution (more advanced version)
//n(r)=ngas(r)*ncr(r),
//n_gas(r)=n1*(1+(r/r0)**2)**(-3*beta/2),
//n_cr(r)=n2*r**(-a).
class HotGasDistributionTimesCR2 : public GalacticDensityDistribution, public DiskCenterCS, protected Parameters
{
public:
	HotGasDistributionTimesCR2(double aSingularitySmoothSizeKpc):
			fMaxRadiusKpc(100),
			fMinRadiusKpc(0),
			fBeta(0.5),
			fR0kpc(3.),
			fR1kpc(0.),//default value 0 left here just for compatibility with previous version
			fA(1.),
			fSingularitySmoothSize(aSingularitySmoothSizeKpc),
			f_D_power(1./3.)//,
			//f_n0(0)
	{
		Reader()->readDoublePar("HG_MaxRadius_kpc", fMaxRadiusKpc);
		Reader()->readDoublePar("HG_MinRadius_kpc", fMinRadiusKpc);
		Reader()->readDoublePar("HG_Beta", fBeta);
		Reader()->readDoublePar("HG_R0_kpc", fR0kpc);
//      In Galaxy mode ProtonsConc parameter now has meaning of protons concentration at Sun's position
//		So we don't edit it anymore
//		f_n0 = 1.35*pow(fR0kpc,-3.*fBeta)/units.cm3;// from arxiv/1412.3116v2
//		IParWriter* w = IParWriter::Instance();
//		if(w)
//			w->writePar("ProtonsConc",f_n0*n_gas(DiskCenterCS::R_sun)*units.cm3);//save the overwritten value in cm^-3
		double pConcSun = 0;
		Reader()->readDoublePar("ProtonsConc", pConcSun);
		pConcSun /= units.cm3;
		//f_n0 =
		std::string customGasDistr = "";
		Reader()->readStringPar("HG_CustomGasDistr", customGasDistr);
		if(customGasDistr.length()>0){
			customGasDistr = "galaxy/" + customGasDistr;
			fCustomGasDistr = new CLinearFunc(new CTableReader(customGasDistr, 2));
			CopyFileOrDir(customGasDistr, plt_local_c);
		}
		bool useCustomCRdistribution = false;
		Reader()->readBoolPar("HG_useCustomCRdistr", useCustomCRdistribution);
		if(useCustomCRdistribution){
			std::string customCRdistr = "";
			Reader()->readStringPar("HG_CustomCRdistr", customCRdistr);
			if(customCRdistr.length()==0){
				fCustomCRdistr = new DefaultCustomDistr(this);
			}
			else {
				std::string tableFile = "galaxy/" + customCRdistr;
				fCustomCRdistr = new MatrixFunction(tableFile);
				CopyFileOrDir(tableFile, plt_local_c);
			}
			CParticleList::Instance()->SetExternal(EProton);//tmp (this is done for calculation of secondaries from pp in halo) todo: read from settings
		}
		else{
			Reader()->readDoublePar("HG_R1_kpc", fR1kpc);
			Reader()->readDoublePar("HG_a", fA);
			double D_cm2_s=0.;//1e29;
			Reader()->readDoublePar("HG_D", D_cm2_s);
			D_atGeV_in_kpc = D_cm2_s*units.cm*units.cm/units.sec/units.Mpc*1e3;
			Reader()->readDoublePar("HG_D_power", f_D_power);
			if(D_atGeV_in_kpc>0){
				//calculate HG_MaxRadius_kpc based on
				double Rmax = sqrt(GalaxyAgeKpc*D_atGeV_in_kpc*pow(Ranges().Emax()/units.GeV, f_D_power));
				if(Rmax<fMaxRadiusKpc){
					std::cerr << "HG_MaxRadius_kpc adjusted to " << Rmax << "using diffusion limits" << std::endl;
					fMaxRadiusKpc = Rmax;
					if(IParWriter::Instance())
						IParWriter::Instance()->writePar("HG_MaxRadius_kpc", fMaxRadiusKpc);
				}
			}
		}
	}

	double MaxDistanceKpc(double l, double b) const
	{
		return GetDistanceFromSun(l, b, fMaxRadiusKpc);
	}

	bool IsIsotropic() const{
		return true;
	}

	double Density(double l, double b, double rKpc) const
	{
		if(fCustomCRdistr)
			return 0.;//GetExternalSpectra() is used
		double r = GetDistanceFromCenter(l, b, rKpc);
		if(r<fMinRadiusKpc)
			return 0.;
		if(r<fSingularitySmoothSize)
			r = fSingularitySmoothSize;
		double result = n_gas(r) * n_cr(r);
		ASSERT_VALID_NO(result);
		return result;
	}

	virtual void GetExternalSpectra(Concentrations& aN, double l, double b, double rKpc, const SourceTable* aSource) const
	{
		if(fCustomCRdistr){
			double r = GetDistanceFromCenter(l, b, rKpc);
			WeightCustom wc(fCustomCRdistr, r);
			//to properly take into account p-p interactions we multiply proton concentration by n_gas(r)/n_gas(sun) and
			//use constant proton density n_gas(sun) for jacobian calculation
			//This procedure also provides correct normalization of output proton spectrum at Sun's position
			aSource->MomentSource(aN, 0., true, &wc, n_gas(rKpc)/n_gas(DiskCenterCS::R_sun)*SourceTable::DefaultMomentSourceDt);
		}
	}

	virtual IFunction* EnergyWeightHotGasDistributionTimesCR(double l, double b, double rKpc) const {
		if(fCustomCRdistr || D_atGeV_in_kpc<=0 || f_D_power==0.)
			return 0;
		double r = GetDistanceFromCenter(l, b, rKpc);
		double Dfrac=r*r/GalaxyAgeKpc/D_atGeV_in_kpc;//Rmax ~ sqrt (t*D)
		double Emin = pow(Dfrac, 1./f_D_power)*units.GeV;//assuming D~E^(1/3)
		return new WeightE(Emin);
	}

	/*
	Medium* CreateMedium() const {
		Medium* m = new Medium();
		m->protonsConc = f_n0*n_gas(DiskCenterCS::R_sun);//set proton concentration to the value corresponding to Sun's position
		return m;
	}*/
private:
	class WeightE : public IFunction
	{
	public:
		WeightE(double aEmin):fEmin(aEmin){};
		virtual double f(double aE) const{
			return aE>fEmin?1.:0.;
		}
	private:
		double fEmin;
	};
	class DefaultCustomDistr : public IFunction2{
	public:
		DefaultCustomDistr(HotGasDistributionTimesCR2* aParent):fParent(aParent){}
		virtual double f(double E, double R) const{
			return fParent->n_cr(R);
		}
	private:
		HotGasDistributionTimesCR2* fParent;
	};

	class WeightCustom : public IFunction
	{
	public:
		WeightCustom(const IFunction2*	aCustomDistr, double aRkpc):
				fCustomDistr(aCustomDistr),fR(aRkpc){
		}
		virtual double f(double aE) const{
			return fCustomDistr->f(aE/units.eV, fR);
		}
	private:
		const IFunction2*	fCustomDistr;
		double fR;
	};

protected:
	virtual double n_cr(double r) const
	{
		double r_smooth = (fR1kpc>0)?(1.+r/fR1kpc):r;
		return pow(r_smooth, -fA);
	}

	virtual double n_gas(double r) const
	{//(1+(r/r0)**2)**(-3*beta/2) here we don't multiply by constant f_n0 (only r-dependence is important)
		double result = 0.;
		if(fCustomGasDistr.isNull()){
			r /= fR0kpc;
			result = pow(1.+r*r,-1.5*fBeta);// from arxiv/1412.3116v2 eq. (1)
		}
		else{
			result = fCustomGasDistr->f(r);
		}
		ASSERT_VALID_NO(result);
		return result;
	}
	double fMaxRadiusKpc;
	double fMinRadiusKpc;
	double fBeta;
	double fR0kpc;//gas distribution plato size
	double fR1kpc;//CR distribution plato size
	double fA;
	double fSingularitySmoothSize;
	double D_atGeV_in_kpc;//diffusion coeff for E=1 GeV in kpc
	double f_D_power;
	//double f_n0;//gas density norm
	static const double GalaxyAgeKpc;//Age of galaxy in kpc
	SafePtr<IFunction2>	fCustomCRdistr;
	SafePtr<IFunction>	fCustomGasDistr;
};

const double HotGasDistributionTimesCR2::GalaxyAgeKpc=3e6;//Age of galaxy in kpc

GalaxyPropagEngine::GalaxyPropagEngine(CInjectionSpectra* aInjection, double aMicroStep):
fInjection(aInjection),
fAngleStep(1.),
fVerboseOutput(false)
{
	redshift.setZ(0.);
	//run before calcCoef since it may modify some media parameters used in calcCoef
	fDensity = GalacticDensityDistribution::Create(Reader(), aMicroStep/units.kpc*0.5);
	fPropagCoef.medium = fDensity->CreateMedium();
	fPropagCoef = calcCoef(fPropagCoef, aMicroStep);
	if(useCELredshift())
		fPropagCoef.jacobian->AddRedshift(0.);
	fSource = new SourceTable(fInjection);

	Reader()->readDoublePar("GalAngleStep", fAngleStep);
    Reader()->readBoolPar("GalVerboseOutput", fVerboseOutput);
}
/*
GalaxyPropagEngine::GalaxyPropagEngine(GalaxyPropagEngine& aOrig):
fOwnsPropagCoef(false),
fInjection(aOrig.fInjection),
fAngleStep(aOrig.fAngleStep),
fVerboseOutput(aOrig.fVerboseOutput)
{
	fPropagCoef.jacobian = aOrig.fPropagCoef.jacobian;
	fPropagCoef.medium = aOrig.fPropagCoef.medium;
	fSource = new SourceTable(fInjection);
	fDensity = initDensity();
}*/

GalaxyPropagEngine::~GalaxyPropagEngine()
{
	fPropagCoef.free();
}

GalacticDensityDistribution* GalacticDensityDistribution::Create(IParReader* aReader, double aSmoothSizeKpc){
	GalDistribution distr = Stellar;
	aReader->readSwitchPar("GalSourceDistribution",&distr,EndGalDistribution);
	switch(distr)
	{
		case Stellar:
			return new StellarDistribution();
		case NFW:
			return new NFWDistribution(aSmoothSizeKpc);
		case HotGas:
			return new HotGasDistributionTimesCR(aSmoothSizeKpc);
		case HotGas2:
			return new HotGasDistributionTimesCR2(aSmoothSizeKpc);
        case GDConst:
            return new ConstGalDistribution();
		case GDBurkert:
			return new BurkertDistribution();
		default:
			ThrowError("Unsupported value of GalSourceDistribution parameter");
	}
	return 0;
}

void GalaxyPropagEngine::run(Concentrations& aN1, double l, double b)//calculate flux from direction l, b (in Galactic coordinates)
{
	double step = fPropagCoef.jacobian->TimeScale();
	double distance = fDensity->MaxDistanceKpc(l,b)*units.Mpc_cm /units.Lunit*1e-3;
	if(step>=distance*0.1)
		ThrowError("GalaxyPropagEngine::run() error: invalid step size");
	fCurL = l;
	fCurB = b;
#ifdef _DEBUG
	cerr << "distance: " << distance * units.Lunit / units.Mpc_cm * 1e3 << "kpc" << endl;
#endif
	CPropagEngine::run(aN1, distance, step, fPropagCoef.jacobian);
}

void GalaxyPropagEngine::run(Concentrations& aN1, double lMin, double lMax, double bMin, double bMax)//calculate average flux
{
	ASSERT(bMax>bMin);
	ASSERT(lMax>lMin);
	//const double omegaTot = (lMax-lMin)*M_PI/180.*fabs(sin(bMin*M_PI/180.)-sin(bMax*M_PI/180.));
	size_t stepsB = round((bMax-bMin)/fAngleStep);
	if(stepsB<1)
		stepsB = 1;
	const double deltaB = (bMax-bMin)/stepsB;
	size_t stepsL = round((lMax-lMin)/fAngleStep);
	if(stepsL<1)
		stepsL = 1;
	const double deltaL = (lMax-lMin)/stepsL;
	double lCur = lMin + 0.5*deltaL;
	double bCur = bMin + 0.5*deltaB;
	Concentrations result;

	Concentrations N;
	double l,b;
	while(true)
	{

		if(lCur>lMax)
		{
			bCur += deltaB;
			lCur = lMin + 0.5*deltaL;
		}
		l = lCur;
		b = bCur;
		lCur += deltaL;

		if(b>bMax)
			break;
		double deltaOmega = deltaL*M_PI/180.*fabs(sin((b+0.5*deltaB)*M_PI/180.)-sin((b-0.5*deltaB)*M_PI/180.));
		N.Copy(aN1);
		run(N, l, b);
        if(fVerboseOutput){
            std::string dir = "galaxy_l" + ToString(l) + "_b" + ToString(b);
            N.PrintSpectrum(dir.c_str());
        }
		N.Multiply(deltaOmega);

		result.Add(N);

	}

	aN1.Copy(result);
}

bool GalaxyPropagEngine::beforeStep(int stepNo, double t, double dt)
{
	double Rkpc = (-0.5*dt-t)*units.Lunit/units.Mpc_cm *1e3;
	double w = fDensity->Density(fCurL, fCurB, Rkpc);
	fCurWeightE = fDensity->EnergyWeight(fCurL, fCurB, Rkpc);
	fSource->SetWeight(w, fCurWeightE);
	fDensity->GetExternalSpectra(*fN,fCurL, fCurB, Rkpc, fSource);
	if(fDensity->UpdateMedium(*fPropagCoef.medium,fCurL, fCurB, Rkpc))
		fPropagCoef.jacobian->Recalculate(*fPropagCoef.medium);
#ifdef _DEBUG
	//cerr << fCurB << "\t" << fCurL << "\t" << Rkpc << "\t" << w << endl;
#endif
	return true;
}

bool GalaxyPropagEngine::start()
{
	fSource->calculateQ();
	return true;
}

GalaxyEffectiveSourceEngine::GalaxyEffectiveSourceEngine(CInjectionSpectra* aInjection, double aMicroStep):
		fInjection(aInjection)
{
	redshift.setZ(0.);
	fPropagCoef = calcCoef(fPropagCoef, aMicroStep);
	if(useCELredshift())
		fPropagCoef.jacobian->AddRedshift(0.);
	fSource = new SourceTable(fInjection);
	fDensity = GalacticDensityDistribution::Create(Reader(), 0.5*aMicroStep/units.kpc);
	if(!fDensity->IsIsotropic())
		ThrowError("GalaxyEffectiveSourceEngine: anisotropic distributions are not supported yet");
}


GalaxyEffectiveSourceEngine::~GalaxyEffectiveSourceEngine()
{
	fPropagCoef.free();
}

void GalaxyEffectiveSourceEngine::run(Concentrations& aN1)
{
	//we assume that density depends on distance from galactic center only
	//and calculate flux using trajectory from galaxy center towards direction opposite to sun (i.e. b=l=0)
	fR = (fDensity->MaxDistanceKpc(0,0)-DiskCenterCS::R_sun)*units.Mpc_cm /units.Lunit*1e-3;

	double step = fPropagCoef.jacobian->TimeScale();
	if(step>=fR*0.1)
		ThrowError("GalaxyEffectiveSourceEngine::run() error: invalid step size");
	CPropagEngine::run(aN1, fR, step, fPropagCoef.jacobian);

	aN1.Multiply(1.2613453127e-05);// normalize for use as custom injection spectrum with CustomInjectionDeltaT_Mpc=1 assuming 1 galaxy per Mpc^-3
    //the norm was obtained by running constant source without interactions
}


bool GalaxyEffectiveSourceEngine::beforeStep(int stepNo, double t, double dt)
{
	//we assume that density depends on distance from galactic center only
	//and calculate flux using trajectory from galaxy center towards direction opposite to sun (i.e. b=l=0)
	double rFromCenterKpc = (fR + t - 0.5*dt)/units.kpc;
	double rFromSunKpc = DiskCenterCS::R_sun + rFromCenterKpc;

	double w = rFromCenterKpc*rFromCenterKpc*fDensity->Density(0, 0, rFromSunKpc);

	fCurWeightE = fDensity->EnergyWeight(0, 0, rFromSunKpc);
	fSource->SetWeight(w, fCurWeightE);
#ifdef _DEBUG
	//cerr << fCurB << "\t" << fCurL << "\t" << Rkpc << "\t" << w << endl;
#endif
	return true;
}

bool GalaxyEffectiveSourceEngine::start()
{
	fSource->calculateQ();
	return true;
}