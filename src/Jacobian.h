#pragma once

#include "Concentrations.h"
#include <gsl/gsl_blas.h>
#include "Vector.h"
#include "Units.h"
#include "Parameters.h"
#include <string>

/// This Jacobian class was primarily designed to be used for rectlinear propagation mode
/// The diffusion approximation is implemented by multiplying all particle rates by (D(E/Z)/D(E_0))^{-1}
/// for rigidity E/Z<E_0 and not modified for E/Z>E_0, where E_0 is free parameter and the diffusion coefficient
/// D(E/z) dependence on rigidity E/Z may have power law form D~(E/Z)^beta or any custom function.
/// Rectlinear propagation mode (beta=0) is used by default
class Jacobian : public GslODEbinningAdaptor,  Parameters
{
public:
	Jacobian(double aTimescale = 1.);
	void SetDiffusion(double aDiffusionBeta, double aDiffusionEmax){
		SetDiffusion((aDiffusionBeta==0) ? 0 : new DefaultD(aDiffusionBeta), aDiffusionEmax);
	}
	void SetLeakyBoxDiffusion(double aDiffusionBeta, TParticle aDiffusionNormParticle, double aFracTesc_Tint, double aNormEnergy, double aLeakyBoxEffPropagTime){
		SetLeakyBoxDiffusion((new DefaultD(aDiffusionBeta)), aDiffusionNormParticle, aFracTesc_Tint, aNormEnergy, aLeakyBoxEffPropagTime);
	}
	//Takes ownership of first argument
	void SetDiffusion(IFunction* aD, double aDiffusionEmax){
		fLeakyBoxDiffusion = false;
		fDiffusion = aD;
		fDiffusionNorm = 1.;
		if(aD){
			fDiffusionNorm = fDiffusion->f(aDiffusionEmax);
			fDiffusionEmax = aDiffusionEmax;
		}
		if(fDiffusionEmax<fMidE[0])
			fDiffusion = 0;
	}

	void SetLeakyBoxDiffusion(IFunction* aD, TParticle aDiffusionNormParticle, double aFracTesc_Tint, double aNormEnergy, double aLeakyBoxEffPropagTime){
		//the particular choice of propag time (aLeakyBoxEffPropagTime param) should not effect final result since only t_esc/t_int metters
		fLeakyBoxDiffusion = true;
		fDiffusionNormParticle = aDiffusionNormParticle;
		fDiffusion = aD;
		fDiffusionNorm = 1.;
		fDNormEnergyBin = Ranges().midE().getBin(aNormEnergy/ParticleData::GetEnergyScaleFactor(aDiffusionNormParticle));
		fFracTesc_Tint=aFracTesc_Tint;
		fLeakyBoxEffPropagTime = aLeakyBoxEffPropagTime;
	}

	double CalculateMaximalRate(TParticle aParticle=EEndParticle) const;

	void Recalculate(const Medium& aCoef);
	void Reset();
	void InterpolateLinearly(const Jacobian& aJ1, const Jacobian& aJ2, double aPart);
	void AddRedshift(double aZ);
	void Create();
	void SetGslJacobian(double *aDfDy) const;
	void PrintInteractionRates(std::string aFileName) const;
	void GetInteractionRate(TParticle aParticle, CVector& aRate) const;
	void PrintSlice(std::string aFileName, TParticle aPrimary, TParticle aSecondary) const;
	//CEL
	static void AddCEL(CMatrixAddOnlyView& coef, const CVector& aEnergyLossRate, int iMinInteractionBin);
	static void AddCELbin(CMatrixAddOnlyView& aCoef, const CVector& aEnergyLossRate, int iPrim);
	static double GetCELintRate(const CVector& aEnergyLossRate, int iPrim);
	
	
	//uses old technic (Implisit scheme)
	void FastEvolve(Concentrations& aConc, const CMatrix* aSource, double aSourceWeight, double aDeltaT) const;
	//todo: move source from PropagCoef
	void SetDirivatives(const double y[], double f[], const CMatrix* aSource, double aSourceWeight) const;
	//calculate maximal source component (multiplied by E) used for source norm ajustment
	static double GetMaxSource(const CMatrix& aSource);
	inline double TimeScale() const {return fTimescale;}
	inline bool RedshiftAdded() const { return fRedshiftAdded;}
    double GetMaxInteractionRate(int max_ebin) const;
private:
	class DefaultD : public IFunction{
	public:
		DefaultD(double beta):fBeta(beta){};
		virtual double f(double E_z)const{ return pow(E_z,fBeta);}
	private:
		double fBeta;
	};
	double fDiffusionNorm;
	double fDiffusionEmax;
	TParticle fDiffusionNormParticle;
	SafePtr<IFunction> fDiffusion;
	bool 	fLeakyBoxDiffusion;
	static int CEL_scheme;
	void AddConstCEL(double aEnergyLossRate, int iMinInteractionBin);

	double												fTimescale;
	bool												fRedshiftAdded;
	CAutoDeletePtrArray<CAutoDeletePtrArray<CMatrix> >	fBonds;
	CVector												fFastEvolveBuffer;
	double	fFracTesc_Tint;
	int 	fDNormEnergyBin;
	double	fLeakyBoxEffPropagTime;
	bool fLimitToInteractionRange;
};
