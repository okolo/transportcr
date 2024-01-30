#ifndef PROPAG_ENGINE_TIPA_ALREADY_INCLUDED
#define PROPAG_ENGINE_TIPA_ALREADY_INCLUDED

#include "Concentrations.h"
#include "Source.h"
#include <gsl/gsl_odeiv.h>
#include "Coupling.h"

class Medium;
class Jacobian;
class Concentrations;

struct TPropagCoef
{
	inline TPropagCoef():medium(0),jacobian(0){}
	inline TPropagCoef(Medium*	aPropagCoef, Jacobian*	aJacobian):medium(aPropagCoef),jacobian(aJacobian){}
	void free();
	Medium*	medium;
	Jacobian*	jacobian;
};

struct ODEParam
{
	const SourceTable*		source;
	const Jacobian*		jacobian;
	double				sourceWeight;
};

enum SourceType
{
	ContinuousSource = 0,
	PointSource
};

class CPropagEngine : protected Parameters
{
public:
	enum TOdeStepMethod
	{//test result on ICS-PP for ODESolvingAccuracy = 1e-4 L=1Mpc microstep=0.1 (made on HP 6730b in VS 2008 debugging mode)
		EFastStep = 0,//old method (1st order implisit scheme with Richardson extrapolation (but without step adjustment)
		EGsl_rk2,//2 sec small problem, bad for ICS point source
		EGsl_rk4,//4 sec ok, bad for ICS point source
		EGsl_rkf45,//2 sec small problem, bad for ICS point source
		EGsl_rkck,//3 sec small problem, bad for ICS point source
		EGsl_rk8pd,//4 sec ok, bad for ICS point source
		EGsl_rk2imp,//6 sec small problem, bad for ICS point source
		EGsl_rk2simp,//uses jacobian, ok for ICS point source (on raw ics 0.1 Mpc 415 sec)
		EGsl_rk4imp,//7 sec ok, bad for ICS point source
		EGsl_bsimp,//uses jacobian (longer than EGsl_rk2simp, same results), ok for ICS point source
		EGsl_gear1,//11 sec ok, bad for ICS point source
		EGsl_gear2,//7 sec ok, bad for ICS point source

		EEndOdeStepMethod
	};

	CPropagEngine();
	virtual ~CPropagEngine();
	void run(Concentrations& aN1, double aDistance, double aMicroStep, Jacobian* aJac1, Jacobian* aJac2=0);
	static inline bool useCELredshift() { return fCELredshift; }
	static SourceType GetSourceType();
protected:
	static void ReadSettings();
	void close();// may be used to enforce reinitialization

	// called after each step
	virtual bool onEndStep(int stepNo, double t, const Jacobian& aJacobianUsed){return true;/*continue calculation*/};
	
	// called before each step
	virtual bool beforeStep(int stepNo, double t, double dt){return true;/*continue calculation*/};
	
	// called before starting calculating cycle
	virtual bool start(){return true;/*continue calculation*/};
	
	// called after finishing calculating cycle
	virtual void onFinish();

	virtual	TPropagCoef calcCoef(TPropagCoef c, double aTimescale);
	
	inline double	getDistance() const{return m_distance;};
	inline double	getMicroStep() const{return dt;};
	inline int		getNoOfSteps() const{return noOfSteps;};
	
	SourceTable*			fSource;
	Concentrations* fN;
private:
	void step(Concentrations& aConc, double _dt, ODEParam& aParam);
    void step_impl(Concentrations& conc, double _dt, ODEParam& aParam);
	void fastStep(double _dt, Concentrations& aConc, ODEParam& aParam);
	void init(int aDim); // lazy initialization
	static const gsl_odeiv_step_type* getStepper();
private:
	double		dt;
	int			noOfSteps;
	double		m_distance;
	static		TOdeStepMethod fOdeStepMethod;
	static double fAbsODESolvingError;
	static double fRelODESolvingError;
	static bool	fCELredshift;
	static bool fReadSettings;

	gsl_odeiv_step*		fGsl_odeiv_step;
	gsl_odeiv_control*	fGsl_odeiv_control;
	gsl_odeiv_evolve*	fGsl_odeiv_evolve;
	//CVector				fFastStepBuffer[2];
	Concentrations* fFastStepBuffer[3];
protected:
	static SourceType		fSourceType;
	bool					fCoefTestEnabled;
	static bool				fEnableDebugOutput;
    static bool             fRuntime_step_correction;
    static double           fStepCorrectionQuantile;
};

class CLongDistanceEngine : protected CPropagEngine
{
public:
	/**
	* Constructor takes ownership of aZbinning
	* aZbinning[0] = 0, aZbinning[n] > aZbinning[n-1]
	*/
	CLongDistanceEngine(const CVector& aZbinning, CInjectionSpectra* aInjection);
	virtual ~CLongDistanceEngine();
	TPropagCoef runZ(Concentrations& aN, double _microStep);
	void setZstep(int _tst);
protected:
	virtual bool isNeeded(TPropagCoef c){return false;};
	virtual bool beforeZstep(){return true;/*continue cycle*/};
	virtual void onZFinish();

	// For use inside beforeStep()
	// Note that here _t is negative
	inline double AgeOfUniverse(double _t){return iGlobalTime + getDistance() + _t;}
/// Overridden from base class
	virtual bool beforeStep(int stepNo, double t, double dt);
	virtual void onFinish();
	virtual bool onEndStep(int stepNo, double t, const Jacobian& aJacobianUsed);

protected:
	double iGlobalTime;
	Concentrations* mp_N;
	///binning on z: mZbinning[0] = 0, mZbinning[n] > mZbinning[n-1]
	const CVector&		mZbinning;
	inline int CurrentZbin() {return mCurrentZbin;}
private:
	int				mCurrentZbin;
	//double			iPrevZ;
	double			iPrevL;
	//do not extrapolate on z propagation coeffitients during calculation of propagation from large z
	static bool NO_Z_COEF_EXTRAPOLATION;

	SmartPtr<CInjectionSpectra> fInjection;
};


class CShortDistanceTestEngine : public CPropagEngine, CouplingParameters
{
public:
	CShortDistanceTestEngine(CInjectionSpectra* aInjection);
	void run(Concentrations& aN1, double aDistance, double aMicroStep);
	void run(Concentrations& aN1, double aMicroStep);
	static double DefaultDistance();
protected:
	virtual bool onEndStep(int stepNo, double t, const Jacobian& aJacobianUsed);
	virtual bool beforeStep(int stepNo, double t, double dt);
	virtual bool start();
	virtual void onFinish();
//private:
	int iOutput;//iterator used for evaluating output time, using output period
	int iMaxOutput;//output period
	int curOutputNo;//current output number
private:
	double iPrevZ;
	TPropagCoef fPropagCoef;
	static double	fDistanceMpc;//distance in Mpc (used if fTauDistance>0 and fTauParticle>=0)
	static int		fNoutputs;//number of intermediate outputs
	static double	fTauDistance;//distance in units of attenuation length
	static int		fTauParticle;//particle used for attenuation length unit of fTauDistance parameter
	static double	fTauEnergyMeV;//Energy in MeV at which the attenuation length is calculated (if<0) the lowest attenuation length is used
	static double	fDiffusionBeta;//diffusion coef power law dependence on energy
	static double   fDiffusionEmax_eV;//maximal energy E/z of diffusion
	static double	fNeutrinoClustering;
	static double	lNeutrinoClustering;
	static bool		fLeakyBoxTest;
	static bool		fSettingsRead;
	SmartPtr<CInjectionSpectra> fInjection;
};

#endif //#ifndef PROPAG_ENGINE_TIPA_ALREADY_INCLUDED
