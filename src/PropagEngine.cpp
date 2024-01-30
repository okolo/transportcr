#include "PropagEngine.h"
#include "Units.h"
#include "TimeZ.h"
#include <math.h>
#include "InjectionSpectra.h"
#include "main.h"
#include "FilePtr.h"
#include "ParticleList.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "Concentrations.h"
#include "Jacobian.h"
#include <time.h>
#include "Coupling.h"
#include "Parameters.h"
#include "Medium.h"
#include <iomanip>

#define SHORT_TEST_DIR "fin"
#define SHORT_TEST_DIR_INI "step0"

int derivatives (double t, const double y[], double f[], void *params)
{
	ODEParam* par = (ODEParam*)params;
//#pragma message("Initialize currentSource")
	par->jacobian->SetDirivatives(y, f, par->source, par->sourceWeight);//if currentSource==0 the source is switched off
	return GSL_SUCCESS;
}

//int callCounter = 0;
int jacobian (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{//check was done if there is need to assign the same values each time
 //the answer is yes
	//callCounter++;
	//cerr << callCounter << endl;
	ODEParam* par = (ODEParam*)params;
	par->jacobian->SetGslJacobian(dfdy);
	gsl_vector v = gsl_vector_view_array(dfdt, par->jacobian->Dim()).vector;
	gsl_vector_set_zero(&v);
	return GSL_SUCCESS;
}

void TPropagCoef::free()
{
	delete jacobian;
	delete medium;
	jacobian = 0;
	medium = 0;
}

CPropagEngine::CPropagEngine():
fSource(0)
{
	CouplingParameters couplingParams;
	fCoefTestEnabled = couplingParams.CoeffTestEnabled();
	fGsl_odeiv_step = 0;
	fGsl_odeiv_control = 0;
	fGsl_odeiv_evolve = 0;
	for(int i=0; i<3; i++)
		fFastStepBuffer[i]=0;

	if(!fReadSettings)
	{
		ReadSettings();
	}
}

void CPropagEngine::ReadSettings()
{
	Reader()->readSwitchPar("OdeStepMethod",&fOdeStepMethod, EEndOdeStepMethod);
	Reader()->readSwitchPar("MomentSource", &fSourceType, 2);
	Reader()->readDoublePar("AbsODESolvingError",fAbsODESolvingError);
	Reader()->readDoublePar("RelODESolvingError",fRelODESolvingError);
	Reader()->readBoolPar("CELredshift",fCELredshift);
	Reader()->readBoolPar("DebugOutput",fEnableDebugOutput);
    Reader()->readBoolPar("Runtime_step_correction", fRuntime_step_correction);
    Reader()->readDoublePar("StepCorrectionQuantile",fStepCorrectionQuantile);
	fReadSettings = true;
}

SourceType CPropagEngine::GetSourceType()
{
	if(!fReadSettings)
	{
		ReadSettings();
	}
	return fSourceType;
}

CPropagEngine::~CPropagEngine()
{
	close();
}

void CPropagEngine::close()
{
	if(fGsl_odeiv_evolve)
		gsl_odeiv_evolve_free (fGsl_odeiv_evolve);
	if(fGsl_odeiv_control)
		gsl_odeiv_control_free (fGsl_odeiv_control);
	if(fGsl_odeiv_step)
		gsl_odeiv_step_free (fGsl_odeiv_step);
	fGsl_odeiv_step = 0;
	fGsl_odeiv_control = 0;
	fGsl_odeiv_evolve = 0;
	delete fFastStepBuffer[1];
	delete fFastStepBuffer[2];
	fFastStepBuffer[1] = 0;
	fFastStepBuffer[2] = 0;
	fFastStepBuffer[0] = 0;
}

double CPropagEngine::fAbsODESolvingError = 1e-4;
double CPropagEngine::fRelODESolvingError = 1e-4;
bool CPropagEngine::fCELredshift = true;
bool CPropagEngine::fEnableDebugOutput = false;
bool CPropagEngine::fRuntime_step_correction = false;
double CPropagEngine::fStepCorrectionQuantile = 0.999;
bool CPropagEngine::fReadSettings = false;
SourceType CPropagEngine::fSourceType = ContinuousSource;

const gsl_odeiv_step_type* CPropagEngine::getStepper()
{
	switch(fOdeStepMethod)
	{
		case EGsl_rk2: return gsl_odeiv_step_rk2;
		case EGsl_rk4: return gsl_odeiv_step_rk4;
		case EGsl_rkf45: return gsl_odeiv_step_rkf45;
		case EGsl_rkck: return gsl_odeiv_step_rkck;
		case EGsl_rk8pd: return gsl_odeiv_step_rk8pd;
		case EGsl_rk2imp: return gsl_odeiv_step_rk2imp;
		case EGsl_rk2simp: return gsl_odeiv_step_rk2simp;
		case EGsl_rk4imp: return gsl_odeiv_step_rk4imp;
		case EGsl_bsimp: return gsl_odeiv_step_bsimp;
		case EGsl_gear1: return gsl_odeiv_step_gear1;
		case EGsl_gear2: return gsl_odeiv_step_gear2;
		default:
			ASSERT(false);
			return 0;
	}
}

// lazy initialization
void CPropagEngine::init(int aDim)
{
	if(fGsl_odeiv_step)
		return;

	fGsl_odeiv_step = gsl_odeiv_step_alloc (getStepper(), aDim);
	fGsl_odeiv_control = gsl_odeiv_control_yp_new (fAbsODESolvingError, fRelODESolvingError);
	fGsl_odeiv_evolve = gsl_odeiv_evolve_alloc (aDim);
}

CPropagEngine::TOdeStepMethod CPropagEngine::fOdeStepMethod = CPropagEngine::EFastStep;

void CPropagEngine::onFinish()
{

}

void CPropagEngine::run(Concentrations& aN1, double aDistance, double aMicroStep, Jacobian* aJac1, Jacobian* aJac2)
{
	//PropagCoef tempCoef;
	fN = &aN1;
	Jacobian interpolatedJacobian;
	bool interpolate = (aJac2 && aJac1!=aJac2);

	m_distance = aDistance;
		
	double t=-m_distance;
	noOfSteps = Round(aDistance/aMicroStep);
	dt=-t/noOfSteps;

	if(!start())
		return;

	ODEParam param;
	param.source = fSource;
	//#pragma message (__FILE__ "(" STRING(__LINE__) "): remove test code")
	param.jacobian = interpolate ? &interpolatedJacobian : aJac1;//permanent code
	//param.jacobian = interpolate ? aJac2 : aJac1;//test code
	double absoluteTime = redshift.t();
	for(int currentStep = 0; currentStep < noOfSteps; currentStep++)
	{
		if(fEnableDebugOutput)
			cerr << "step " << currentStep << " of " << noOfSteps << endl;

		if(!beforeStep(currentStep,t,dt))
			break;

		//Concentrations conc(*fN);

		//before solving dif equation adjust concentration norm to make max(source*dt, concentration) to be of order of unity
		double concNorm = aN1.GetMaxValue();
		double sourceNorm = fSource?Jacobian::GetMaxSource(*fSource)*dt:0;
		double norm = concNorm > sourceNorm ? concNorm : sourceNorm;

		if(fEnableDebugOutput)
			cerr << "z=" << redshift.z() << "\t\tconc=" << concNorm << "\t\tsource=" << sourceNorm << endl;

		ASSERT(norm>=0);
		if(norm>0)//propagate only if concentration is non-zero or/and source is non-zero
		{
			if(interpolate)
			{
				interpolatedJacobian.InterpolateLinearly(*aJac1, *aJac2, ((double)currentStep)/((double)noOfSteps));
				if(useCELredshift())
					interpolatedJacobian.AddRedshift(CTimeZ::t2z(absoluteTime));
			}
			aN1.Multiply(1./norm);
			param.sourceWeight = 1./norm;
			step(aN1, dt, param);
			aN1.Multiply(norm);//recover original norm
			//conc.Copy(*fN);
			//nonJacobianInteractions(dt);
		}
		t+=dt;
		absoluteTime+=dt;
		if(!onEndStep(currentStep,t,*param.jacobian))
			break;
	}//end for
	onFinish();
}

void CPropagEngine::fastStep(double aDt, Concentrations& aConc, ODEParam& aParam)
{
	fFastStepBuffer[0] = &aConc;
	int ac,ntst;
	for(ac=1;ac<3;ac++)
	{
		if(fFastStepBuffer[ac]==0)
			fFastStepBuffer[ac] = new Concentrations();
		ASSERT(fFastStepBuffer[ac-1]->Dim()==aConc.Dim());
		fFastStepBuffer[ac]->Copy(aConc);
	}
		
	for(ac=0,ntst=2;ac<3;ac++,ntst+=2)
	{
		double deltaT=aDt/ntst;
		for(int step=ntst;step>0;step--)
		{
			aParam.jacobian->FastEvolve(*(fFastStepBuffer[ac]), aParam.source, aParam.sourceWeight, deltaT);
		};
	};

	int iMax = aConc.Dim();
	double x0=aDt/2.0;
	double x1=aDt/4.0;
	double x2=aDt/6.0;
	double A0=x1*x2/(x2-x0)/(x1-x0);
	double A1=x0*x2/(x2-x1)/(x0-x1);
	double A2=x1*x0/(x0-x2)/(x1-x2);

	for(int i=0; i<iMax; i++)
	{
		double& conc0 = (fFastStepBuffer[0])->Data()[i];
		const double conc1 = (fFastStepBuffer[1])->Data()[i];
		const double conc2 = (fFastStepBuffer[2])->Data()[i];
		conc0 = A0*conc0 + A1*conc1 + A2*conc2;
		TESTVALUE(conc0, conc2);
	};
	fFastStepBuffer[0] = 0;
}

void CPropagEngine::step_impl(Concentrations& conc, double _dt, ODEParam& aParam)
{
    const double rescaledDt = _dt/aParam.jacobian->TimeScale();
    ASSERT((aParam.jacobian->RedshiftAdded() && useCELredshift()) || ((!aParam.jacobian->RedshiftAdded()) && !useCELredshift()));

    if(fOdeStepMethod == EFastStep)
    {
        fastStep(rescaledDt, conc, aParam);
    }
    else
    {
        init(conc.Dim());// lazy initialization

        gsl_odeiv_system sys = {derivatives, jacobian, (size_t)conc.Dim(), &aParam};

        double t = 0., t1 = rescaledDt;
        double h = 0.5*(t1-t);//1e-6

        while (t < t1)
        {
            int status = gsl_odeiv_evolve_apply (fGsl_odeiv_evolve, fGsl_odeiv_control, fGsl_odeiv_step,
                                                 &sys,
                                                 &t, t1,
                                                 &h, conc.Data());

            if (status != GSL_SUCCESS)
                throw "ODE solving error";
        }
        conc.TrimNegative(fAbsODESolvingError);
    }
}

void CPropagEngine::step(Concentrations& conc, double _dt, ODEParam& aParam)
{
    if (!fRuntime_step_correction){
        step_impl(conc, _dt, aParam);
        return;
    }

    while(_dt > 0){
        int max_substeps = 1000;
        int max_ebin = conc.GetEnergyQuantileBin(fStepCorrectionQuantile);
        double q_energy = Ranges().midE()[max_ebin] / units.eV;
        double max_rate = aParam.jacobian->GetMaxInteractionRate(max_ebin);
        double step = max_rate > 0 ? 1./max_rate : _dt;
        double tot_n_steps = floor(_dt / step + 0.5); // using double since tot_n_steps may be beyond int32 or even int64 range
        if (tot_n_steps < 1.0)
            tot_n_steps = 1.0;
        step = _dt / tot_n_steps;
        int n_steps = max_substeps;
        if (((double)n_steps) > tot_n_steps){
            n_steps = int(tot_n_steps + 0.1);
            _dt = 0; // end cycle after this iteration
        }
        else
            _dt -= (n_steps*step);
        std::cerr << "step scale factor " <<  tot_n_steps;
        std::cerr << " for Emax = " << std::setprecision(2) << std::scientific << q_energy << " eV" << std::endl;

        for(int i=0; i < n_steps; i++)
            step_impl(conc, step, aParam);
    }
}

CLongDistanceEngine::CLongDistanceEngine(const CVector&	aZbinning, CInjectionSpectra* aInjection):
mp_N(NULL),
mZbinning(aZbinning),
mCurrentZbin(-1),
iPrevL(-1),
fInjection(aInjection)
{
	READ_BOOL_SETTING(NO_Z_COEF_EXTRAPOLATION);
	if(fSourceType == ContinuousSource)
		fSource=new SourceTable(fInjection);
#ifdef _DEBUG
	ASSERT(useCELredshift()||((aZbinning.length()-1)%Ranges().AccuracyZ()==0))
	ASSERT(mZbinning[0]==0);
	ASSERT(mZbinning.length()>=2);
	for(int i=1; i<mZbinning.length(); i++)
	{
		ASSERT(mZbinning[i]>mZbinning[i-1]);
	}
#endif

}

CLongDistanceEngine::~CLongDistanceEngine()
{
	delete fSource;
}


void CLongDistanceEngine::setZstep(int _tst)
{
	mCurrentZbin = _tst;
	redshift.setZ(mZbinning[mCurrentZbin]);
}

TPropagCoef CLongDistanceEngine::runZ(Concentrations& aN, double aMicroStep)
// returns PropagCoef at z=0
{
	if(aMicroStep<0)
		ThrowError("Invalid step size");
	mp_N = &aN;
	
	setZstep(mZbinning.length()-1);
	iPrevL = 1. + mZbinning[mCurrentZbin];
	//double Zmax = CInjectionSpectra::getZcalculationMax();

	TPropagCoef coef1, coef2;
	coef1 = calcCoef(coef1, aMicroStep);

	if(fSourceType == PointSource)
	{
		SourceTable s(fInjection);
		s.MomentSource(aN, redshift.z(), false);
		fSource=0;
		aN.PrintSpectrum("iniL");
	}
	else
		fSource = new SourceTable(fInjection);
	for(; mCurrentZbin>0; setZstep(mCurrentZbin-1))
	{
		if(!beforeZstep())
			break;
		
		if(redshift.z()<Ranges().Zmax())
		{
			LOG_MESSAGE4("ntst=",mCurrentZbin," z=", redshift.z());
		};
		double zNext=mZbinning[mCurrentZbin - 1];
		iGlobalTime = redshift.t();
		double deltaT = CTimeZ::z2t(zNext) - iGlobalTime;
		double theMicroStep = aMicroStep/pow(1.+redshift.z(),3);//CMB concentration and therefore proc speed is proportional to (1+z)^3
		{
			int ntst = (int)(deltaT/theMicroStep);
			if (ntst<2) {
				ntst = 2;
			}
			theMicroStep = deltaT/ntst;
			LOG_MESSAGE2("number of microsteps ", ntst);
		}
		
		/// calculate coefficients for the next redshift
		redshift.setZ(mZbinning[mCurrentZbin-1]);	
		coef2 = calcCoef(coef2, aMicroStep);
		redshift.setZ(mZbinning[mCurrentZbin]);
		
		Jacobian* j_final=(NO_Z_COEF_EXTRAPOLATION)?NULL:coef2.jacobian;		

		run( aN, deltaT, theMicroStep, coef1.jacobian, j_final);

		if((!useCELredshift())&&(mCurrentZbin-1)%Ranges().AccuracyZ()==0)
		{
			aN.Expansion();
			iPrevL = 1. + mZbinning[mCurrentZbin - 1];
		}

		if(fCoefTestEnabled){
			char dirName[20];
			sprintf(dirName,"Zstep%d",mCurrentZbin);
			aN.PrintSpectrum(dirName);
		}
			
		TPropagCoef mem = coef1;
		coef1 = coef2;
		coef2 = mem;
	}

	Mkdir(plt_local_dir + "uniform");

	if(!isNeeded(coef2))
	{
		delete coef2.jacobian;
		delete coef2.medium;
	}

	onZFinish();
	return coef1;
}

bool CLongDistanceEngine::NO_Z_COEF_EXTRAPOLATION = false;//do not extrapolate on z propagation coeffitients during calculation of propagation from large z

void CLongDistanceEngine::onZFinish()
{
}

TPropagCoef CPropagEngine::calcCoef(TPropagCoef c, double aTimescale)
{
	if(c.medium==NULL)
	{
		c.medium = new Medium;
	}
	if(c.jacobian==NULL)
		c.jacobian = new Jacobian(aTimescale);
	ASSERT(c.jacobian->TimeScale()==aTimescale);
	
	if(c.medium->Update())
		throw "Coef calculation error";
	c.jacobian->Recalculate(*(c.medium));
	if(fCoefTestEnabled)
	{
		CouplingList::Instance()->DebugOutput(plt_local_dir.c_str());
		c.jacobian->PrintInteractionRates(plt_local_dir + "ratesZ=" + ToString(redshift.z()));
		c.jacobian->PrintSlice(plt_local_dir + "gamma_gamma_Z=" + ToString(redshift.z()), EPhoton, EPhoton);
	}
	return c;
}

bool CLongDistanceEngine::beforeStep(int stepNo, double t, double dt)
{
	double age = AgeOfUniverse(t);
	double z = CTimeZ::t2z(age);
	double energyCorrectionFactor = 1.;

	if(iPrevL>=1. && !CPropagEngine::useCELredshift())
	{
		energyCorrectionFactor = (1.+z)/iPrevL;//<=1 energyCorrectionFactor was introduced to compensate discontinuous way of red shifting (since rev 65)
	}

	if(fSource!=NULL)
		fSource->calculateQ(AgeOfUniverse(t+0.5*dt), energyCorrectionFactor);// note that here t is negative
	if(fCoefTestEnabled)
        fN->AppendRecordToEnergyLog(age * units.Lunit / units.Mpc, ToString(z));//energy test
	return true;
}

double s_outputTravelDistancesMpc[] = {-1.}; //{7.68194, 47.7795, 68.0223, 77.9939, 79.6257, 79.7352, 138.597, 138.977, 139.389, 139.787, 183.649, 235.935, 236.102, -1.}; //elements must be in accending order except the last one, which must be less than zero
static int s_nextOutput = 0;

bool CLongDistanceEngine::onEndStep(int stepNo, double t, const Jacobian& aJacobianUsed){
    if(s_outputTravelDistancesMpc[s_nextOutput] < 0)
        return true;

    double now = AgeOfUniverse(t);
    double start = CTimeZ::z2t(mZbinning[mZbinning.length()-1]);
    double distanceTravelled = (now-start)/units.Mpc;
    if (distanceTravelled>s_outputTravelDistancesMpc[s_nextOutput]){
        double z = CTimeZ::t2z(now);
        fN->PrintSpectrum(ToString(s_outputTravelDistancesMpc[s_nextOutput]) + "z" + ToString(z));
        s_nextOutput++;
    }
    return true;
}

void CLongDistanceEngine::onFinish()
{
	close(); // closing concentrations since next time they will be
	//created with different index shift
}

bool CShortDistanceTestEngine::start()
{
	LOG_MESSAGE("Short distance test:");
	if(fSourceType == ContinuousSource)
	{
		fN->Reset();
		fSource = new SourceTable(fInjection);
		fSource->calculateQ();
	}
	else
	{
		SourceTable s(fInjection);
		s.MomentSource(*fN, 0., false);
		fSource = 0;
		fN->PrintSpectrum(SHORT_TEST_DIR_INI);
		if(fCoefTestEnabled)
            fN->AppendRecordToEnergyLog((-getDistance()) * units.Lunit / units.Mpc);//energy test
	}
	iPrevZ = CTimeZ::d2z(getDistance());
		
	iOutput = 0;//iterator used for evaluating output time, using output period
	iMaxOutput = getNoOfSteps()/fNoutputs;//output period
	curOutputNo = 1;//current output number
	return true;
}

void CShortDistanceTestEngine::run(Concentrations& aN1, double aDistance, double aMicroStep)
{
	if(!fPropagCoef.jacobian)
		fPropagCoef.jacobian = new Jacobian(aMicroStep);
	if(fLeakyBoxTest){
		fPropagCoef.jacobian->SetLeakyBoxDiffusion(fDiffusionBeta,(TParticle)fTauParticle, fTauDistance, fTauEnergyMeV*units.MeV, aDistance);
	}else {
		fPropagCoef.jacobian->SetDiffusion(fDiffusionBeta, fDiffusionEmax_eV * units.eV);
	}
	fPropagCoef = calcCoef(fPropagCoef, aMicroStep);
	time_t t = time(0);
	if(useCELredshift()){
		if(fDiffusionEmax_eV * units.eV > Ranges().Emin() && fDiffusionBeta != 0)
			WARN("Redshift disabled : not supported in diffusion mode");
		else
			fPropagCoef.jacobian->AddRedshift(0.);
	}
	CPropagEngine::run(aN1, aDistance, aMicroStep, fPropagCoef.jacobian);
	t = (time(0)-t);
	LOG_MESSAGE3("Short distance test Calculation took ", t, " sec");
	delete fPropagCoef.jacobian;
	delete fPropagCoef.medium;
	fPropagCoef.jacobian = 0;
	fPropagCoef.medium = 0;
}

CShortDistanceTestEngine::CShortDistanceTestEngine(CInjectionSpectra* aInjection):
iPrevZ(-1.),
fInjection(aInjection)
{
	if(!fSettingsRead)
	{
		READ_DOUBLE_SETTING(fNeutrinoClustering);
		READ_DOUBLE_SETTING(lNeutrinoClustering);
		fSettingsRead = true;
	}
	fPropagCoef.jacobian = 0;
	fPropagCoef.medium = 0;

	Reader()->readDoublePar("L_short_distance_test", fDistanceMpc);
	Reader()->readIntPar("shortDistanceTestOutputNo", fNoutputs);
	Reader()->readDoublePar("sdtTauDistance", fTauDistance);
	Reader()->readIntPar("sdtTauParticle", fTauParticle);
	Reader()->readDoublePar("sdtTauEnergyMeV", fTauEnergyMeV);
	Reader()->readDoublePar("sdtDiffusionBeta", fDiffusionBeta);
	Reader()->readDoublePar("sdtDiffusionEmax", fDiffusionEmax_eV);
	Reader()->readBoolPar("LeakyBoxTest", fLeakyBoxTest);
}

void CShortDistanceTestEngine::run(Concentrations& aN1, double microStep)
{
	double distance = DefaultDistance();
	double step = microStep;
	if(fTauParticle>=0 && fTauDistance>0 && !fLeakyBoxTest)
	{
		bool COEF_TEST_ONsave = fCoefTestEnabled;
		fCoefTestEnabled = false;
		TPropagCoef c;
		c.jacobian = new Jacobian(microStep);
		c.jacobian->SetDiffusion(fDiffusionBeta, fDiffusionEmax_eV * units.eV);
		c = calcCoef(c, microStep);
		fCoefTestEnabled = COEF_TEST_ONsave;
		CVector intRates;
		c.jacobian->GetInteractionRate((TParticle)fTauParticle, intRates);
        double maxRate =c.jacobian->CalculateMaximalRate();
        delete c.jacobian;
		delete c.medium;
		double rate = 0.;
		double tauEnergyMeV = fTauEnergyMeV;
		if(tauEnergyMeV>0.)
		{
			CLinearFunc intRateF(Ranges(), intRates);
			intRateF.SetExtension(ExtConst);
			rate = intRateF(tauEnergyMeV/units.Eunit);
		}
		else
		{//finding maximal interaction rate
			int iSelected = -1;
			for(int i = intRates.length()-1; i>=0; i--)
				if(rate<intRates[i])
				{
					rate = intRates[i];
					iSelected = i;
				}
			if(rate>0)
				tauEnergyMeV = Ranges().midE()[iSelected]*units.Eunit;
		}
		if(rate<=0.)
		{
			ThrowError(
					ToString("Failed to set short distance test length based on interaction length: interaction length is 0 for particle ")
					+ ParticleData::getParticleFileName((TParticle)fTauParticle));
		}
		logger.Write(
				ToString("Using interaction length ") +	ToString(1./rate/units.Mpc_cm *units.Lunit) + " Mpc obtained for " +
				ParticleData::getParticleFileName((TParticle)fTauParticle) + " at energy " + ToString(tauEnergyMeV*1e6) + " eV"
				);
		//setting distance based on interaction rate
		distance = fTauDistance/rate;

		//adjusting step size
        if(step==0.){//automatic stepping adjustment
            step = 1./maxRate;
            logger.Write(ToString("Setting step to minimal interaction length ") +	ToString(step/units.Mpc) + " Mpc");
        }
        if(step<0.)
        {//step was given in units of tau and not Mpc
            double stepInTau = -step/units.Mpc;
            logger.Write(ToString("Adjusting step to ") + ToString(stepInTau)  + " tau = " + ToString(stepInTau/rate/units.Mpc) + " Mpc");
            step = stepInTau/rate;
        }
        else if(step>0.1/rate){
            step = 0.1/rate;
            logger.Write(ToString("Step adjusted to ") + ToString(step/units.Mpc) + " Mpc");
        }
	}
	//adjusting step size for proper intermediate output
	if(step>distance/fNoutputs)
		step = distance/fNoutputs;
	run(aN1, distance, step);
}

double	CShortDistanceTestEngine::fDistanceMpc = 100.;//distance in Mpc (used if fTauDistance>0 and fTauParticle>=0)
int		CShortDistanceTestEngine::fNoutputs = 1;//number of intermediate outputs
double	CShortDistanceTestEngine::fTauDistance = 1;//distance in units of attenuation length
int		CShortDistanceTestEngine::fTauParticle = -1;//particle used for attenuation length unit of fTauDistance parameter
double	CShortDistanceTestEngine::fTauEnergyMeV = -1.;//Energy in MeV at which the attenuation length is calculated (if<0) the lowest attenuation length is used
double CShortDistanceTestEngine::fDiffusionEmax_eV = 0;//Maximal energy E/z for which diffusion mode is used
double CShortDistanceTestEngine::fDiffusionBeta = 0.333333333333333;
double CShortDistanceTestEngine::fNeutrinoClustering=1;//neutrino clustering factor (typical 100-300). Setting to 1 switches off the clustering effect
double CShortDistanceTestEngine::lNeutrinoClustering=0;//here in Mpc, further in internal units (typical 5). Setting to 0 switches off the clustering effect
bool CShortDistanceTestEngine::fLeakyBoxTest = false;
bool CShortDistanceTestEngine::fSettingsRead = false;


double CShortDistanceTestEngine::DefaultDistance()
{
	return fDistanceMpc*units.Mpc_cm /units.Lunit;
}

bool CShortDistanceTestEngine::onEndStep(int stepNo, double t, const Jacobian& aJacobianUsed)
{
	if(fCoefTestEnabled)
	{
		if(stepNo == 0)
		{
			aJacobianUsed.PrintInteractionRates(plt_local_dir + "ratesZ=" + ToString(redshift.z()));
			aJacobianUsed.PrintSlice(plt_local_dir + "gamma_gamma_Z=" + ToString(redshift.z()), EPhoton, EPhoton);
		}
        fN->AppendRecordToEnergyLog(t * units.Lunit / units.Mpc, "actual_E.plt");//energy test
	}

	iOutput++;
	if((iOutput>=iMaxOutput)&&stepNo)
	{
		double distanceTraveled = (getDistance()+t)*units.Lunit/units.Mpc_cm;
		if(distanceTraveled>getDistance())
			distanceTraveled = getDistance();
		char dirName[128];
		
		sprintf(dirName,"step%d",curOutputNo);
		fN->PrintSpectrum(dirName);
		CFilePtr distFile(Fopen(plt_local_dir+dirName+DIR_DELIMITER_STR+"distanceTraveled","wt"));
		fprintf(distFile,"%lg",distanceTraveled);
		iOutput=0;
		curOutputNo++;
	}
	return true;
}

bool CShortDistanceTestEngine::beforeStep(int stepNo, double t, double dt)
{
	//shortDistanceTestTime = t*Lunit/Mpc;

	double prev = Medium::neutrinoClusteringModifier;
	Medium::neutrinoClusteringModifier=(-t<lNeutrinoClustering)?fNeutrinoClustering:1.;
	if(
		(prev != Medium::neutrinoClusteringModifier && ((!NO_NEUTRINO_s) || (!NO_NEUTRINO_t))))
	{
		fPropagCoef.medium->Update();
		fPropagCoef.jacobian->Recalculate(*fPropagCoef.medium);
	}
	LOG_MESSAGE2("step ", stepNo);
	return true;
}

void CShortDistanceTestEngine::onFinish()
{
	LOG_MESSAGE("done");
	fN->PrintSpectrum(SHORT_TEST_DIR);
	fPropagCoef.medium->background()->saveBackground(SHORT_TEST_DIR);
	if(fCoefTestEnabled)
        fN->AppendRecordToEnergyLog(0);//energy test
}
