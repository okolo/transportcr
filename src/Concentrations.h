#pragma once

#include "Vector.h"
#include "TParticle.h"
#include "Ranges.h"

//class Concentrations;

class GslODEbinningAdaptor
{
public:
	inline int Dim() const {return fDim;}
protected:
	GslODEbinningAdaptor();

	//returns negative value if aParticle is not propagating
	inline	int	Index(TParticle aParticle, int iE) const
	{
		return fParticleIndex[aParticle] * fNumberOfEnergyBins + iE;
	}

	inline	int	RawIndex(int aPropagParticle, int iE) const
	{
		return aPropagParticle * fNumberOfEnergyBins + iE;
	}
	const CBinning&			fMidE;
	int						fNumberOfEnergyBins;
	const Array<int,int&>&	fParticleIndex;
	int						fNumberOfPropagatingParticles;
	int						fDim;
};

class Concentrations : public GslODEbinningAdaptor
{
public:
	Concentrations();
	Concentrations(const Concentrations& aConcentr);
	Concentrations(double* aData);
	~Concentrations(void);

	void ReadAll(std::string aFilePath);
	double TotalEnergy(TParticle aParticle = EEndAllParticles);
	int AppendRecordToEnergyLog(double aCurT/*current value of evol parameter (typically time or redshift)*/,
                                std::string aExtraData = "",/* added to the end of the record*/
                                std::string aFileName = "energy");// used for debuging
	void Reset();
	void Expansion();
	void Copy(const Concentrations& c);
	void PrintSpectrum(std::string _dir);
	double GetValue(TParticle aParticle, double aE) const;

	//get maximal value of vector component
	double GetMaxValue() const;
	//multiply vector by constant (used for normalization)
	void Multiply(double aFactor);
	//add another vector to this
	void Add(const Concentrations& aConcentr);
	void SetZero();
	bool IsValid();
	void AssertValid();
	void TrimNegative(double aAccuracy);
    int GetEnergyQuantileBin(double quantile) const;

	inline double Bin(TParticle aParticle, int iE) const
	{
		int index = Index(aParticle, iE);//returnes negative value if aParticle is not propagating
		return index < 0 ? 0. : fData[index]/fMidE[iE];
	}
	inline void SetBin(TParticle aParticle, int iE, double aValue)
	{
		//Assertion fails if aParticle is not propagating
		fData[Index(aParticle, iE)] = fMidE[iE] * aValue;
	}
	inline double& Data(int aParticleIndex, int iE)
	{
		return fData[aParticleIndex * fNumberOfEnergyBins + iE];
	}

	inline double Data(int aParticleIndex, int iE) const
	{
		return fData[aParticleIndex * fNumberOfEnergyBins + iE];
	}

	inline double* Data(){return fData;}

	inline const double* Data() const {return fData;}

private:
	double*					fData;
	bool					fOwnsData;
};

