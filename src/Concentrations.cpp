#include "Concentrations.h"
#include "ParticleList.h"
#include "gsl/gsl_vector_double.h"
#include "FilePtr.h"
#include "TableReader.h"
#include "TableFunc.h"
#include "ParticleData.h"
#include "Units.h"
#include "Coupling.h"

GslODEbinningAdaptor::GslODEbinningAdaptor():
fMidE(Ranges().midE()),
fParticleIndex(CParticleList::Instance()->GetPropagatingParticleIndex()),
fNumberOfPropagatingParticles(CParticleList::Instance()->NumberOfPropagatingParticles())
{
	fNumberOfEnergyBins = Ranges().nE();
	fDim = fNumberOfEnergyBins * fNumberOfPropagatingParticles;
}

Concentrations::Concentrations():
	fData(0),
	fOwnsData(true)
{
	fData = new double[Dim()];
	Reset();
}

Concentrations::Concentrations(const Concentrations& aConcentr):
	fData(0),
	fOwnsData(true)
{
	fData = new double[Dim()];
	Copy(aConcentr);
}

void Concentrations::Copy(const Concentrations &c)
{
	ASSERT(Dim()==c.Dim());
	memcpy(fData, c.fData, sizeof(double)*Dim());
}

void Concentrations::Reset()
{
	memset(fData,0,sizeof(double)*Dim());
}

/*
Concentrations::Concentrations(const Concentrations& aConcentr):
	fData(0),
	fOwnsData(true)
{
	fData = new double[Dim()];
	FOR_ALL_PARTICLES_INVOLVED(particle)
	{
		const CParticle& p = aConcentr.getParticle((TParticle)particle);
		int index = Index(particle, 0);//returnes negative value if aParticle is not propagating
		ASSERT(index >= 0);
		for(int i=0; i<fNumberOfEnergyBins; i++)
			fData[index+i] = p[i]*fMidE[i];
	}
}
*/

void Concentrations::ReadAll(std::string aFilePath)
{
	CFilePtr in(Fopen(aFilePath, "rt"));
	if((FILE*)in == NULL)
		ThrowError("Failed to load spectrum from " + aFilePath);

	CTableReader reader(in,EEndAllParticles + 2);
	CVector x(reader.getColumn(0));
	double coeff_y = 1./( (BC().ss2-BC().ss_2)*units.Vunit*4.0*Pi/units.lightspeed/units.Eunit*units.outEunitMeV);
	double coeff_x = units.Eunit/units.outEunitMeV;
	FOR_ALL_PARTICLES_INVOLVED(particle)
	{
		CVector y(reader.getColumn(particle + 1));
		CDefaultTableFunc func(x,y);
		double scaleE = ParticleData::GetEnergyScaleFactor(particle);
		for(int i=0; i<fNumberOfEnergyBins; i++)
		{
			double E = fMidE[i]*scaleE;
			SetBin(particle, i, func.f(E*coeff_x)/E/coeff_y);
		}
	}
}

int Concentrations::GetEnergyQuantileBin(double quantile) const
{
    CouplingParameters cp;
    CVector energy(fNumberOfEnergyBins);
    FOR_ALL_PARTICLES_INVOLVED(particle){
        if(!cp.Interacts(particle))
            continue; // we count only interacting particles
        int index = Index(particle, 0);
        double sum = 0;
        for(int i = 0; i<fNumberOfEnergyBins; i++){
            sum += fData[index + i]*fMidE[i];
            energy[i] += sum;
        }
    }
    double max_energy = quantile * energy[fNumberOfEnergyBins-1];
    for(int i = 0; i<fNumberOfEnergyBins; i++)
        if (energy[i] > max_energy)
            return i;

    return fNumberOfEnergyBins - 1;
}

void Concentrations::Expansion()
{
	FOR_ALL_PARTICLES_INVOLVED(particle)
	{
		for(int i=1;i<fNumberOfEnergyBins;i++)
			SetBin((TParticle)particle,i-1,Bin((TParticle)particle,i));
		SetBin((TParticle)particle,fNumberOfEnergyBins-1,0);
	}
}

double Concentrations::TotalEnergy(TParticle aParticle)
{
	double energy=0;
    if(aParticle==EEndAllParticles)//calculate total energy in all particles
    {
        for (int i = 0; i < fNumberOfEnergyBins; i++) {
            double sum = 0.;
            FOR_ALL_REAL_PARTICLES_INVOLVED(particle) {
                if (particle == EStartNuclei)
                    break;
                sum += Bin((TParticle) particle, i);
            }
            FOR_ALL_NUCLEI_INVOLVED(nucleus) {
                double multiplier = nucleus - EStartNuclei + 2;
                double flux = Bin((TParticle) nucleus, i);
                sum += (multiplier * flux);
            }
            energy += fMidE[i] * sum;
        }
    }
    else{//calculate energy in just one type of particle
        double scalingFactor = ParticleData::GetEnergyScaleFactor(aParticle);
        for(int i=0; i<fNumberOfEnergyBins; i++)
        {
            energy += fMidE[i]*Bin(aParticle, i)*scalingFactor;
        }
    }
	return energy;
}

int Concentrations::AppendRecordToEnergyLog(double x, std::string aExtraData, std::string fileName)
{
	CFilePtr ff(FopenL(fileName, plt_local_c,"a"));
	if(ff==0)
	{
		cerr << "failed to open " << fileName;
		return -1;
	}
    fprintf(ff,"\n%lg",x);
    double totalE = 0.;
    FOR_ALL_PARTICLES(particle)//print values even for particles being excluded from calculation to preserve file format
    {
        const char* particleName = ParticleData::getParticleFileName((TParticle)particle);
        if(particleName){
            double energy = TotalEnergy((TParticle)particle)*units.Eunit/units.outEunitMeV/units.Vunit;
            fprintf(ff,"\t%lg",energy);
            totalE += energy;
        }
    };
    fprintf(ff,"\t%lg",totalE);
    if(aExtraData.length())
        fprintf(ff,"\t%s", aExtraData.c_str());
	return 0;
}

void Concentrations::SetZero()
{
	gsl_vector v = gsl_vector_view_array(fData, Dim()).vector;
	gsl_vector_set_zero(&v);
}

Concentrations::Concentrations(double* aData):
fData(aData),
fOwnsData(false)
{
	
}

Concentrations::~Concentrations(void)
{
	if(fOwnsData)
		delete[] fData;
	fData = 0;
}

/*
void Concentrations::Copy(class Concentrations & aConcentr) const
{
	FOR_ALL_PARTICLES_INVOLVED(particle)
	{
		CParticle& p = aConcentr.getParticle((TParticle)particle);
		int index = Index(particle, 0);//returnes negative value if aParticle is not propagating
		ASSERT(index >= 0);
		for(int i=0; i<fNumberOfEnergyBins; i++)
			p[i] = fData[index+i]/fMidE[i];
	}
}*/

double Concentrations::GetValue(TParticle aParticle, double aE) const
{
	double binE = aE/ParticleData::GetEnergyScaleFactor(aParticle);
	int iE = Ranges().midE().findLeftX(binE);
	if(iE<0)
		return 0;
	double leftVal = Bin(aParticle,iE);
	if(iE==fNumberOfEnergyBins-1)
		return leftVal;
	double rightVal = Bin(aParticle,iE+1);
	double leftE = Ranges().midE()[iE];
	double deltaE = Ranges().midE()[iE+1] - leftE;
	return leftVal + (rightVal-leftVal)/deltaE*(binE-leftE);
}

void Concentrations::PrintSpectrum(std::string _dir)
{
	const double Emin = Ranges().Emin();
	const double Emax = Ranges().Emax();

	FILE* output=NULL;
	string full_dir = (plt_local_dir + _dir);
	int i;
	double E;
	double coeff_y=1./((BC().ss2-BC().ss_2)*units.Vunit*4.0*Pi/units.lightspeed/units.Eunit*units.outEunitMeV);
	double coeff_x=1./(units.outEunitMeV /units.Eunit);

	LOG_MESSAGE3("Writing data to ",full_dir,"...")

	Mkdir(full_dir);

	output=Fopen("all",full_dir,"wt",1);
	fprintf(output,"#Energy\tj(E)*E^2:");
	//const char* particleNames[EEndParticle];

	{
		int particle;
		for(particle = 0, i=2;particle<(int)EEndAllParticles;i++,particle++)
		{
			const char* particleName = ParticleData::getParticleFileName((TParticle)particle);
			if(particleName)
				fprintf(output,"%d=%s\t",i,particleName);
			else
				i--;
		}
		fprintf(output,"%d=All",i);
	}
	double max = CParticleList::Instance()->NucleiEnabled()?Emax*56:Emax;
	for(i=0, E=Emin*BC().ss2; E<max; E*=BC().ss1, i++)
	{
		double sumAll = 0.;
		fprintf(output,"\n%lg",E*coeff_x);
		FOR_ALL_PARTICLES(particle)//print values even for particles being excluded from calculation to preserve file format
		{
			const char* particleName = ParticleData::getParticleFileName((TParticle)particle);
			if(particleName){
				double flux = (particle<EStartNuclei && i < fNumberOfEnergyBins) ? Bin((TParticle)particle, i) : GetValue(
						(TParticle) particle, E);
				flux *= E*coeff_y;
				fprintf(output,"\t%lg",flux);
				sumAll += flux;
			}
		};
		fprintf(output,"\t%lg",sumAll);
	};
	fclose(output);

	LOG_MESSAGE("Done");
}

//get maximal value of vector component
double Concentrations::GetMaxValue() const
{
	double result = 0;
	for(int i=Dim()-1; i>=0; i--)
		if(fData[i]>result)
			result = fData[i];
	return result;
}

bool Concentrations::IsValid()
{
	for(int i=Dim()-1; i>=0; i--)
		if(fData[i]<0)
			return false;
	return true;
}

void Concentrations::AssertValid()
{
	for(int i=Dim()-1; i>=0; i--)
		if(fData[i]<0)
		{
			int index = i/fNumberOfEnergyBins;
			TParticle part = CParticleList::Instance()->GetPropagatingParticles()[index];
			string text = "invalid ";
			text += ParticleData::getParticleFileName(part);
			text += "[" + ToString(i-index*fNumberOfEnergyBins) + "]=" + ToString(fData[i]);
			logger.Write(text, Log::EError);
			throw "invalid vector";
		}
}

void Concentrations::TrimNegative(double aAccuracy)
{
	for(int i=Dim()-1; i>=0; i--)
		if(fData[i]<0)
		{
			if(-fData[i] > aAccuracy)
			{
				int index = i/fNumberOfEnergyBins;
				TParticle part = CParticleList::Instance()->GetPropagatingParticles()[index];
				string text = "invalid ";
				text += ParticleData::getParticleFileName(part);
				text += "[";
				text += ToString(i%fNumberOfEnergyBins);
				text += "] flux (negative value exceeding accuracy level " + ToString(aAccuracy) + ")";
				logger.Write(text, Log::EError);
				ASSERT(0);
			}
			fData[i] = 0;
		}
}

//multiply vector by constant (used for normalization)
void Concentrations::Multiply(double aFactor)
{
	for(int i=Dim()-1; i>=0; i--)
		fData[i] *= aFactor;
}


void Concentrations::Add(const Concentrations& aConcentr)
{
	ASSERT(Dim()==aConcentr.Dim());
	for(int i=Dim()-1; i>=0; i--)
		fData[i] += (aConcentr.fData[i]);
}
