/*
 * AA.h
 *
 *  Created on: May 5, 2010
 *      Author: kalashev
 */

#ifndef AA_H_
#define AA_H_
#include "configFortran.h"
#include "Coupling.h"
#include "TableFunc.h"

namespace couplings {

	enum PPCrossSectionMode
{
	PPAnalKelner, //p-p cross sections implementation based on parameterizations of Kelner et al. astro-ph/0606058v1
	PPCodeKachelriess, //M. Kachelriess cross section tables built using QGSJet along with low energy part parametrization by T. Kamae et al., ApJ 647 (2006) 692
	PPTableKachelriess // Using cross section tables only from Kachelriess (still under development)
};

//p-p interaction with energy dependent grammage support
class NucleonProton: public Coupling, public IFunction {
	class Channel_N_x : public CouplingChannelT<NucleonProton>
	{
	public:
		Channel_N_x(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary) : CouplingChannelT<NucleonProton>(aCoupling, aPrimary, aSecondary){};
		void Coef(CMatrixAddOnlyView& aCoef) const;
	protected:
		virtual double F(double x, double L) const = 0;
	};
	class Channel_N_xKachelriess : public CouplingChannelT<NucleonProton>
	{
	public:
		Channel_N_xKachelriess(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary);
		void Coef(CMatrixAddOnlyView& aCoef) const;

		//d(sigma)/dE * E in internal sigma units
		double sigma(double Eprim, double Esec) const;
	private:
		static bool fInitialized;
	};
	class Channel_N_xTable : public CouplingChannelT<NucleonProton>
	{
	public:
		Channel_N_xTable(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary);
		void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		MatrixFunction fSigma;
	};
	class Channel_N_mnu : public Channel_N_x
	{
	public:
		Channel_N_mnu(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary) : Channel_N_x(aCoupling, aPrimary, aSecondary){};
		double F(double x, double L) const {return F_e(x,L)+F_mnu1(x,L);}
	};
	class Channel_N_e_ne : public Channel_N_x
	{
		public:
		Channel_N_e_ne(NucleonProton* aCoupling, TParticle aPrimary, TParticle aSecondary) : Channel_N_x(aCoupling, aPrimary, aSecondary){};
		double F(double x, double L) const {return F_e(x,L); }
	};
	class Channel_N_gamma : public Channel_N_x
	{
	public:
		Channel_N_gamma(NucleonProton* aCoupling, TParticle aPrimary) : Channel_N_x(aCoupling, aPrimary, EPhoton){};
		double F(double x, double L) const { return F_gamma(x,L); }
	};
	class Channel_N_N : public CouplingChannelT<NucleonProton>
	{
	public:
		Channel_N_N(NucleonProton* aCoupling, TParticle aPrimary) : CouplingChannelT<NucleonProton>(aCoupling, aPrimary, aPrimary){};
		void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	virtual double f(double _x) const{
		return Sigma_inel(_x);
	}
	double Sigma_inel(double aEp) const;
	double Sigma_inel_anal(double aEp) const;
	double Sigma_inel_tab(double aEp) const;
	// aTable parameter can be used  to define rates either by optical depth dependence on energy (eV)
	// or by dependence of grammage in g/cm^2 on energy in eV
	// In both cases grammage/depth is converted
	// to effective rates for constant travel time t=1 Mpc for all energies
	// This is useful for calculation of effective spectra of compact sources
	// For this one should run close source mode with L_short_distance_test=1
	NucleonProton(RateMode mode=RMdefault, PPCrossSectionMode aScMode=PPCodeKachelriess, const char* aTable = 0);
	void InitWithTableRate();
	void InitWithTableGrammage();
	static double F_gamma(double x, double L);
	static double F_mnu1(double x, double L);
	static double F_e(double x, double L);
	void SetBackgrounds(const Medium& aPropagCoef);
private:
	double 		fPconc;
	double 		fEth;
	CVector 	fRates;
	std::string	fTableFile;
	RateMode 	fMode;
	PPCrossSectionMode fScMode;
	SafePtr<CLogScaleLinearFunc> fTableSigma;
	double		fTableSigmaLeftMult;
	double		fTableSigmaRightMult;
	static const double		fTableEmin_eV[];
	static const double		fTableEmax_eV[];
};

	// Nuclei with proton gas interaction with energy dependent grammage support
	// Secondaries calculation is not supported yet
	class NucleiProton: public Coupling {
		class Channel_A_A : public CouplingChannelT<NucleiProton>
		{
		public:
			Channel_A_A(NucleiProton* aCoupling, TParticle aPrimary) : CouplingChannelT<NucleiProton>(aCoupling, aPrimary, aPrimary){};
			void Coef(CMatrixAddOnlyView& aCoef) const;
		};
	public:
		//aTable parameter can be used to define dependence of grammage (in g/cm^2) on E/Z
		//this is useful for calculation of effective spectra of compact sources
		//In this case one should run close source mode with L_short_distance_test=1 Mpc
		//(since grammage is converted to effective rates derived for constant travel time t=1Mpc for all energies)
		NucleiProton(const char* aGrammageTable = 0);
		void SetBackgrounds(const Medium& aPropagCoef);
	private:
		CMatrix fRatesPerUnitDensity;
		double fLastConc;
		double fRateMult;
	};

}

#endif /* AA_H_ */
