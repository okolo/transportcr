/*
 * BLLac.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: ok
 */

#include "BLLac.h"
#include "TimeZ.h"
#include "Units.h"

const double BLLacLDDE::Lmin = 1e-10;
const double BLLacLDDE::Lmax = 100.;

BLLacLDDE::BLLacLDDE() : fNL0(0) {
	fNL0 = NL(0.);
}

double BLLacLDDE::e(double z, double L/*in units of 10^48 erg/s*/)
{
	const double p1 = 1.54;
	const double p2 = -0.42;
	const double Zs = 2.1;
	const double alpha = 0.052;

	double Zc = Zs*pow(L,alpha);
	double frac = (1.+z)/(1+Zc);
	double result = 1./(pow(frac,p1) + pow(frac,p2));
	ASSERT_VALID_NO(result);
	return result;
}

double BLLacLDDE::rho(double L/*in units of 10^48 erg/s*/)
{
	const double A_Log10 = 41.7;
	const double Lc = 1.82;
	const double g1 = 0.59;
	const double g2 = 1.43;
	L /= Lc;
	double result = A_Log10/(L*(pow(L,g1) + pow(L,g2)));
	ASSERT_VALID_NO(result);
	return result;
}

void BLLacLDDE::UnitTest()
{//compare with dN/dz from Fig. 2 of arXiv/1311.5708v3
	double H_in_km_s_Mpc = 71;
	double Lv = 0.73;
	CRanges::SetDefault(
			new CRanges(
				10,
				1e6*units.eV,
				1e21*units.eV,
				100,
				1,
				false,
				1.0,
				0
				)
			);
	CTimeZ::init(Lv,H_in_km_s_Mpc);

	double step = pow(10.,0.005);
	double Mpc = units.Mpc_cm /units.Lunit;//Mpc in internal units
	std::cout << "# z\t n(z)\t n_L(z)\tdV_c/dz" << std::endl;
	for(double z1=1; z1<9; z1*=step)
	{
		double z = z1-1;
		double Dc = CTimeZ::z2Dc(z);
		double dV_dz = 4*M_PI*Dc*Dc*CTimeZ::diffDc(z)/(Mpc*Mpc*Mpc);
		//todo: take into account \omega(F_\gamma) - Fermi-LAT efficiency at flux (F_\gamma)
		std::cout << z << "\t" << N(z) << "\t" << NL(z) << "\t" << dV_dz << std::endl;
	}
}

double BLLacLDDE::N(double z)
{
	class NKern : public IFunction
	{
		double fZ;
	public:
		NKern(double z) : fZ(z){}
		double f(double logL) const {
			double L = exp(logL);
			return L*e(fZ,L)*rho(L);
		}
	};

	NKern f(z);
	double result = funcUtils.gsl_integration_qag(&f,log(Lmin),log(Lmax), 0, 1e-4, 10000);
	ASSERT_VALID_NO(result);
	return result;
}

double BLLacLDDE::NL(double z)
{
	class NLKern : public IFunction
	{
		double fZ;
	public:
		NLKern(double z) : fZ(z){}
		double f(double logL) const
		{
			double L = exp(logL);
			return e(fZ,L)*rho(L)*L*L;
		}
	};
	NLKern f(z);
	double result = funcUtils.gsl_integration_qag(&f,log(Lmin),log(Lmax), 0, 1e-4, 10000);
	ASSERT_VALID_NO(result);
	return result;
}
