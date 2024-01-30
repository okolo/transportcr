#include "ParticleData.h"
#include "Medium.h"
#include "Jacobian.h"
#include "Coupling.h"
#include "TimeZ.h"
#include "FilePtr.h"
#include "Units.h"

Jacobian::Jacobian(double aTimescale):
    fTimescale(aTimescale),
    fRedshiftAdded(false),
    fLeakyBoxDiffusion(false),
    fLimitToInteractionRange(false)
{
    Reader()->readBoolPar("LimitToInteractionRange", fLimitToInteractionRange);
}

void Jacobian::SetGslJacobian(double *aDfDyData) const
{//todo: optimize with use of memcpy
	gsl_matrix_view view = gsl_matrix_view_array(aDfDyData,Dim(),Dim());
	gsl_matrix*		DfDyMatrix = &view.matrix;
	gsl_matrix_set_zero(DfDyMatrix);

	const Array<TParticle, TParticle&>& particles = CParticleList::Instance()->GetPropagatingParticles();
	CParticleList* pl = CParticleList::Instance();

	for(int iSec=0; iSec<fNumberOfPropagatingParticles; iSec++)
	{
		TParticle secondary = particles[iSec];
		if(pl->IsExternal(secondary))
			continue;//don't evolve externally given spectra
		const int indexSec = RawIndex(iSec, 0);
		for(int iPrim=0; iPrim<fNumberOfPropagatingParticles; iPrim++)
		{
			CMatrix* bond = fBonds(iSec)[iPrim];
			if(!bond)
				continue;
			const int indexPrim = RawIndex(iPrim, 0);
			for(int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
			{
				for(int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++)
					gsl_matrix_set (DfDyMatrix, indexSec+eSec, indexPrim+ePrim, (*bond)[eSec][ePrim]);
			}
		}
	}
}

void Jacobian::SetDirivatives(const double y[], double f[], const CMatrix* aSource, double aSourceWeight) const
{
	const Concentrations conc((double*)y);
	Concentrations deriv(f);
	const Array<TParticle, TParticle&>& particles = CParticleList::Instance()->GetPropagatingParticles();
	CParticleList* pl = CParticleList::Instance();

	for(int iSec=0; iSec<fNumberOfPropagatingParticles; iSec++)
	{
		TParticle secondary = particles[iSec];
		if(pl->IsExternal(secondary)) {//don't evolve externally given spectra
			for(int eSec = 0; eSec < fNumberOfEnergyBins; eSec++) {
				deriv.Data(iSec, eSec) = 0.;
			}
			continue;
		}

		double weight = aSourceWeight*fTimescale;
		for(int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
		{
			double& der = deriv.Data(iSec, eSec);
			der = aSource ? weight*fMidE[eSec]*(*aSource)[secondary][eSec] : 0.;
			for(int iPrim=0; iPrim<fNumberOfPropagatingParticles; iPrim++)
			{
				CMatrix* bond = fBonds(iSec)[iPrim];
				if(!bond)
					continue;
				
				for(int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++)
					der += (*bond)[eSec][ePrim] * conc.Data(iPrim, ePrim);
			}
		}
	}
}

void Jacobian::FastEvolve(Concentrations& aConc, const CMatrix* aSource, double aSourceWeight, double aDeltaT) const
{
	CParticleList* pl = CParticleList::Instance();
	CVector& buf = (CVector&)fFastEvolveBuffer;
	if(!buf.length())
		buf.create(Dim());
	buf.copy(aConc.Data());
	Concentrations orig(buf.ptr());
	const Array<TParticle, TParticle&>& particles = CParticleList::Instance()->GetPropagatingParticles();

	for(int iSec=0; iSec<fNumberOfPropagatingParticles; iSec++)
	{
		TParticle secondary = particles[iSec];
		if(pl->IsExternal(secondary))
			continue;//don't evolve externally defined spectra
		double weight = aSourceWeight*fTimescale;
        int eSec = fLimitToInteractionRange ? Ranges().nMinInteraction(secondary) : 0;
        for(; eSec < fNumberOfEnergyBins; eSec++)
		{
			double income = aSource ? weight*fMidE[eSec]*(*aSource)[secondary][eSec] : 0.;
			double absorption = fBonds(iSec)[iSec] ? -(*(fBonds(iSec)[iSec]))[eSec][eSec] : 0.;
			for(int iPrim=0; iPrim<fNumberOfPropagatingParticles; iPrim++)
			{
				CMatrix* bond = fBonds(iSec)[iPrim];
				if(!bond)
					continue;

                int ePrim = (iSec == iPrim) ? eSec : Ranges().nMinInteraction(particles[iPrim]);
				for(; ePrim < fNumberOfEnergyBins; ePrim++)
					income += (*bond)[eSec][ePrim] * orig.Data(iPrim, ePrim);
			}
			double origValue = orig.Data(iSec, eSec);
			if(absorption>0)
				income += absorption*origValue;//subtract absorption term from income

			double val = origValue+aDeltaT*income;
			if(val<0)
				throw "Fast evolve scheme failed: income<0";

			VERIFY_VALID_NO(aConc.Data(iSec, eSec) = val/(1.0+aDeltaT*absorption));
		}
	}
}

double Jacobian::GetMaxSource(const CMatrix& aSource)
{
	const CBinning& E = Ranges().midE();
	double result = 0.;
	int nn = Ranges().nE();
	FOR_ALL_REAL_PARTICLES_INVOLVED(particle)
	{
		for(int i=0; i<nn; i++)
		{
			double source = aSource[particle][i]*E[i];
			if(source > result)
				result = source;
		}
	}
	return result;
}

void Jacobian::Create()
{
	ASSERT(!fBonds.length());
	for(int i=0; i<fNumberOfPropagatingParticles; i++)
	{
		fBonds.add(new CAutoDeletePtrArray<CMatrix>());
		fBonds(i).create(fNumberOfPropagatingParticles);
	}
}

void Jacobian::AddRedshift(double aZ)
{
	ASSERT(!fRedshiftAdded);
	ASSERT(fBonds.length() == fNumberOfPropagatingParticles);//Create() or Recalculate(const PropagCoef& aCoef) should be called prior
	fRedshiftAdded = true;
	AddConstCEL(CTimeZ::eLossRate(aZ), 0);
}

enum TCelScheme
{
	CelSchemeStable=0,//standard CEL approximation way: less accurate but more stable (accuracy checked on redshift only mode)
	CelSchemeBeta, //more accurate scheme which may not work for Fast ODE solving scheme
//both modes seem to give the same result in the limit of infinite number of energy bins
	CelSchemeEOF //should be the last (used for schemes number count only)
};

int Jacobian::CEL_scheme = CelSchemeStable;

void Jacobian::AddCEL(CMatrixAddOnlyView& aCoef, const CVector& aEnergyLossRate, int iMinInteractionBin)
{
	const int nn = Ranges().nE();

	for(int iPrim = iMinInteractionBin; iPrim<nn; iPrim++)
		AddCELbin(aCoef, aEnergyLossRate, iPrim);
}

double Jacobian::GetCELintRate(const CVector& aEnergyLossRate, int iPrim)
{
	double multCEL = 1./BC().sLn;	
	if(CEL_scheme == CelSchemeStable)
		return multCEL * aEnergyLossRate[2*iPrim];
	else
		return -multCEL*(aEnergyLossRate[2*iPrim+2] - aEnergyLossRate[2*iPrim] - 1.5*aEnergyLossRate[2*iPrim+1]);
}

double Jacobian::GetMaxInteractionRate(int aMaxBinE) const
{
    if(fBonds.length()==0)
        return 0.;
    double max_rate = 0.;
    FOR_ALL_REAL_PARTICLES_INVOLVED(particle) {
        int index = fParticleIndex[particle];
        if(index < 0)
            continue;
        CMatrix* bond = fBonds(index)[index];
        if(bond)
            for(int i=0; i<=aMaxBinE; i++) {
                double rate = (-(*bond)[i][i]);
                if (rate > max_rate)
                    max_rate = rate;
            }
    }
    return max_rate/fTimescale;
}

void Jacobian::AddCELbin(CMatrixAddOnlyView& aCoef, const CVector& aEnergyLossRate, int iPrim)
{
	double multCEL = 1./BC().sLn;

	if(CEL_scheme == CelSchemeStable)
	{
		double coef = multCEL * aEnergyLossRate[2*iPrim];
		
		aCoef.Add(iPrim, iPrim, -coef);
		if(iPrim>0)
		{
			aCoef.Add(iPrim-1, iPrim, coef);
		}
	}
	else
	{
		double coef0 = (aEnergyLossRate[2*iPrim+2] - aEnergyLossRate[2*iPrim] - 1.5*aEnergyLossRate[2*iPrim+1])*multCEL;
		aCoef.Add(iPrim, iPrim, coef0);
		if(iPrim>0)
		{
			double coef1 = 2.*multCEL*aEnergyLossRate[2*iPrim-1];
			aCoef.Add(iPrim-1, iPrim, coef1);
		}
		if(iPrim>1)
		{
			double coef2 = -0.5*multCEL*aEnergyLossRate[2*iPrim-3];
			aCoef.Add(iPrim-2, iPrim, coef2);
		}
	}
}

void Jacobian::AddConstCEL(double aEnergyLossRate, int iMinInteractionBin)
{
	const int nn = Ranges().nE();
	double multCEL = 1./BC().sLn;

	double coef = multCEL*aEnergyLossRate*fTimescale;
	double coef0 = -1.5*coef;
	double coef1 = 2.*coef*BC().ss_1;
	double coef2 = -0.5*coef*BC().ss_1*BC().ss_1;

	//double eLossLength = 1/(multCEL*aEnergyLossRate)*Lunit/Mpc;//debug

	for(int particle=0; particle<fNumberOfPropagatingParticles; particle++)
	{
		CMatrix* bond = fBonds(particle)[particle];
		if(!bond)
		{
			bond = new CMatrix(fNumberOfEnergyBins);
			fBonds(particle)[particle] = bond;
		}
		for(int iPrim = iMinInteractionBin; iPrim<nn; iPrim++)
		{
			if(CEL_scheme == CelSchemeStable)
			{
				(*bond)[iPrim][iPrim] -= coef;
				if(iPrim>0)
				{
					(*bond)[iPrim-1][iPrim] += (coef*BC().ss_1);
				}
			}
			else
			{
				(*bond)[iPrim][iPrim] += coef0;
				if(iPrim>0)
				{
					(*bond)[iPrim-1][iPrim] += coef1;
				}
				if(iPrim>1)
				{
					(*bond)[iPrim-2][iPrim] += coef2;
				}
			}
		}
	}
}

void Jacobian::Reset()
{
	fRedshiftAdded = false;
	if(!fBonds.length())
		Create();	
	else
		for(int i=0; i<fNumberOfPropagatingParticles; i++)
		{
			for(int j=0; j<fNumberOfPropagatingParticles; j++)
			{
				CMatrix* bond = fBonds(i)[j];
				if(bond)
					bond->reset();
			}
		}
}

void Jacobian::Recalculate(const Medium& aCoef)
{
	Reset();
	CPointerArray<Coupling>& couplings = CouplingList::Instance()->Couplings();
	int nCoupl = couplings.length();
	for(int iCoupl = 0; iCoupl < nCoupl; iCoupl++)
	{
		Coupling& coupl = couplings(iCoupl);
		if(!coupl.IsEnabled())
			continue;
		coupl.SetBackgrounds(aCoef);
		int nChannels = coupl.Channels().length();
		for(int iChannel = 0; iChannel <  nChannels; iChannel++)
		{
			const CouplingChannel& cc = coupl.Channels()(iChannel);
			if(!cc.IsEnabled())
				continue;

			int i = fParticleIndex[cc.Secondary()];
			int j = fParticleIndex[cc.Primary()];
			CMatrix* bond = fBonds(i)[j];
			if(!bond)
			{
				bond = new CMatrix(fNumberOfEnergyBins);
				fBonds(i)[j] = bond;
			}
			CMatrixAddOnlyView mView(*bond);
			cc.Coef(mView);
		}
	}
	///transforming N_i jacobian to (N_i E_i) jacobian with time unit aTimescale
	for(int iSec=0; iSec<fNumberOfPropagatingParticles; iSec++)
		for(int iPrim=0; iPrim<fNumberOfPropagatingParticles; iPrim++)
		{
			CMatrix* bond = fBonds(iSec)[iPrim];
			if(!bond)
				continue;
			for(int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
				for(int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++)
					(*bond)[eSec][ePrim] *= (fTimescale*fMidE[eSec]/fMidE[ePrim]);
		}

	if(fDiffusion.isNull())
		return;

    PrintInteractionRates(plt_local_dir + "rates");
    /// effective rate modification with diffusion
	if(fLeakyBoxDiffusion) { // use formula (C5) of arXiv/1505.02153v2
		int iNormParticle = fParticleIndex[fDiffusionNormParticle];
		double normRate = -((*(fBonds(iNormParticle)[iNormParticle]))[fDNormEnergyBin][fDNormEnergyBin])/fTimescale;
		double normTau = log(1.+fFracTesc_Tint);
		double rateMultNorm = normTau/(normRate*fLeakyBoxEffPropagTime);
		double t_esc_norm = fFracTesc_Tint/(rateMultNorm*normRate);
		double regidity_norm = Ranges().midE()[fDNormEnergyBin]*ParticleData::GetEnergyScaleFactor(fDiffusionNormParticle)/ParticleData::getElectricCharge(fDiffusionNormParticle);
		double diffusionNorm = fDiffusion->f(regidity_norm);

		FOR_ALL_PARTICLES_INVOLVED(primary) {
			int z = abs(ParticleData::getElectricCharge(primary));
			if (!z)
				continue;
			int iPrim = fParticleIndex[primary];
			CMatrix& totRates = *(fBonds(iPrim)[iPrim]);
			double factor = ParticleData::GetEnergyScaleFactor(primary) / z;
			for (int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++) {
				double E_z = fMidE[ePrim] * factor;

				double t_esc = t_esc_norm * diffusionNorm / fDiffusion->f(E_z);
				double orig_rate = -(totRates[ePrim][ePrim])/fTimescale;
                ASSERT_VALID_NO(orig_rate);
                if(orig_rate>0) {
                    double eff_tau = log(1. + t_esc * orig_rate * rateMultNorm);
                    ASSERT_VALID_NO(eff_tau);
                    double mult = eff_tau / (orig_rate * fLeakyBoxEffPropagTime);
                    ASSERT_VALID_NO(mult);

                    for (int iSec = 0; iSec < fNumberOfPropagatingParticles; iSec++) {
                        CMatrix *bond = fBonds(iSec)[iPrim];
                        if (!bond)
                            continue;

                        for (int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
                            (*bond)[eSec][ePrim] *= mult;
                    }
                }
			}
		}

	} else
	{
		FOR_ALL_PARTICLES_INVOLVED(primary) {
			int z = abs(ParticleData::getElectricCharge(primary));
			if (!z)
				continue;
			int iPrim = fParticleIndex[primary];
			CMatrix *totRates = fBonds(iPrim)[iPrim];
			double factor = ParticleData::GetEnergyScaleFactor(primary) / z;
			for (int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++) {
				double E_z = fMidE[ePrim] * factor;
				if (E_z > fDiffusionEmax)
					break;

				double mult = (fDiffusionNorm / fDiffusion->f(
						E_z));/// old model assuming that optical depth is proportional to escape time devided by interaction time

				for (int iSec = 0; iSec < fNumberOfPropagatingParticles; iSec++) {
					CMatrix *bond = fBonds(iSec)[iPrim];
					if (!bond)
						continue;

					for (int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
						(*bond)[eSec][ePrim] *= mult;
				}
			}
		}
	}
    PrintInteractionRates(plt_local_dir + "rates-w-diffusion");
}

void Jacobian::InterpolateLinearly(const Jacobian& aJ1, const Jacobian& aJ2, double aPart)
{
	Reset();
	ASSERT(aJ1.fTimescale == aJ2.fTimescale);
	ASSERT(aJ1.RedshiftAdded() == aJ2.RedshiftAdded());
	fTimescale = aJ1.fTimescale;
	fRedshiftAdded = aJ1.RedshiftAdded();
	if(!fBonds.length())
		Create();

	ASSERT(aPart >= 0 && aPart <= 1);

	double coef1 = 1.-aPart;
	double coef2 = aPart;

	for(int i=0; i<fNumberOfPropagatingParticles; i++)
	{
		for(int j=0; j<fNumberOfPropagatingParticles; j++)
		{
			CMatrix* bond1 = aJ1.fBonds(i)[j];
			if(bond1)
			{
				CMatrix* bond2 = aJ2.fBonds(i)[j];
				ASSERT(bond2);
				CMatrix* bond = fBonds(i)[j];
				if(!bond)
				{
					bond = new CMatrix(fNumberOfEnergyBins);
					fBonds(i)[j] = bond;
				}
				for(int eSec = 0; eSec < fNumberOfEnergyBins; eSec++)
					for(int ePrim = 0; ePrim < fNumberOfEnergyBins; ePrim++)
						(*bond)[eSec][ePrim] = coef1*(*bond1)[eSec][ePrim]+coef2*(*bond2)[eSec][ePrim];
			}
		}
	}
}

void Jacobian::PrintInteractionRates(std::string aFileName) const
{
	CFilePtr file(FopenL(aFileName, "wt"));
	fprintf(file, "#interaction rates in Mpc^-1\n#1-st column energy in eV (E/A for nuclei); n-th column rate for (n-2)-th particle (see TParticle enum)\n");
	const CBinning& E = Ranges().midE();
	const int nn = E.length();
	for(int i = 0; i<nn; i++)
	{
		fprintf(file, "%g", E[i]*units.Eunit*1e6);
		FOR_ALL_PARTICLES(particle)
		{
			int index = fParticleIndex[particle];
			CMatrix* bond = (index < 0 || fBonds.length()==0) ? 0 : fBonds(index)[index];
			if(bond)
			{
				double rate = (-(*bond)[i][i])/fTimescale/units.Lunit*units.Mpc_cm;//interaction rate in Mpc^-1
				fprintf(file,"\t%g", rate);
			}
			else
				fprintf(file,"\t0");
		}
		fprintf(file,"\n");
	}
}

double Jacobian::CalculateMaximalRate(TParticle aParticle /* =EEndParticle */) const{
    double maxRate = 0;
    const CBinning& E = Ranges().midE();
    const int nn = E.length();
    for(int i = 0; i<nn; i++)
    {
        FOR_ALL_PARTICLES(particle)
        {
            if(aParticle==particle || aParticle==EEndParticle){
                int index = fParticleIndex[particle];
                const CMatrix* bond = (index < 0 || fBonds.length()==0) ? 0 : fBonds(index)[index];
                if(bond)
                {
                    double rate = (-(*bond)[i][i])/fTimescale;//interaction rate in internal units
                    if (rate>maxRate)
                        maxRate = rate;
                }
            }
        }
    }
    return maxRate;
}

void Jacobian::GetInteractionRate(TParticle aParticle, CVector& aRates) const
{
	int nn = Ranges().nE();
	if(aRates.length()==0)
		aRates.create(nn);
	ASSERT(aRates.length()==nn);
	int index = fParticleIndex[aParticle];
	CMatrix* bond = (index < 0 || fBonds.length()==0) ? 0 : fBonds(index)[index];
	if(!bond)
		aRates.reset();
	else{
		for(int i = 0; i<nn; i++)
		{
			aRates[i] = (-(*bond)[i][i])/fTimescale;
		}
	}
}

void Jacobian::PrintSlice(std::string aFileName, TParticle aPrimary, TParticle aSecondary) const
{
	std::string sliceName = ParticleData::getParticleFileName(aPrimary);
	sliceName = sliceName + "->" + ParticleData::getParticleFileName(aSecondary);

	const CBinning& E = Ranges().midE();
	const int nn = E.length();
	int iPrim = fParticleIndex[aPrimary];
	int iSec = fParticleIndex[aSecondary];
	CMatrix* bond = (iPrim < 0 || iSec<0 || fBonds.length()==0) ? 0 : fBonds(iSec)[iPrim];
	if(!bond)
	{
		LOG_MESSAGE3("skipping jacobian slice output ", sliceName, " (zero matrix)");
		return;
	}

	CFilePtr file(FopenL(aFileName, "wt"));
	fprintf(file, "#%s interaction rates in Mpc^-1\n#row number - secondary energy bin; column number - primary energy bin\n", sliceName.c_str());

	for(int eSec = 0; eSec<nn; eSec++)
	{
		for(int ePrim = 0; ePrim<nn; ePrim++)
		{
			double rate = ((*bond)[eSec][ePrim])*E[ePrim]/E[eSec]/fTimescale/units.Lunit*units.Mpc_cm;//rate in Mpc^-1
			fprintf(file,"%g\t", rate);
		}
		fprintf(file,"\n");
	}
}
