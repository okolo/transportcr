#if !defined(MEDIUM_H__INCLUDED)
#define MEDIUM_H__INCLUDED
#include "Concentrations.h"
#include "MassiveNeutrino.h"
#include "Background.h"
#include "Vector.h"

class CBackgroundIntegral;

class Medium : Parameters
{
public:
	Medium();
	~Medium(void);

	double neutrinoConc;
	double protonsConc;

	//magnetic field in internal units
	inline double MagnField() const
	{
		return fCurrentB;
	}
	int Update();
	static IBackgroundSpectrum* GetEBL();

	inline const SmartPtr<CBackgroundTable>& background() const {return m_pBackground;};
	inline void setBackground(CBackgroundTable* aBackgound) {m_pBackground =  aBackgound;};

	const CBackgroundIntegral& BackgroundIntegral() const {return *iBackgroundIntegral;}

	static double neutrinoClusteringModifier;//current clustering modifier (=1. no clustering)
private:
	static IBackgroundSpectrum* AddIRBebl();
	static void AddRadioEbl();

	SmartPtr<CBackgroundTable>	m_pBackground;
	CBackgroundIntegral*		iBackgroundIntegral;
	double fB0;
	double fCurrentB;
	bool fBincr;//is magnetic field increasing with redshift
	static SmartPtr<CompoundBackground> fEBL;
	double neutrinoConc0;
	bool	fCoefTestEnabled;
};


#endif//end of file
