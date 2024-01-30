/*
	fundamental consants initialization, and F-classes' members definitions
*/
/*****************************************************************************************/
#include "const.h"
#include "ParticleData.h"
#include "Units.h"

//const double Pi=3.1415927;

//class FWeak - weak interaction constants

double FWeak::G=1.16639e-11;//Fermi constant / MeV^(-2)

//const double FWeak::gL_l=-0.76895;//gL = -cos(ThetaW)^2 //old
const double FWeak::gL_l=-0.73105;//gL = -1/2-sin(ThetaW)^2//Dima
const double FWeak::gR_l=0.23105;//gR = sin(ThetaW)^2
const double FWeak::gL_n=0.5;//gL = 1/2
const double FWeak::gR_n=0.;//gR = 0


double FWeak::Mz=91173.0;//MeV
double FWeak::Mw=80220.0;//MeV
double FWeak::GammaZ=2487.0;//MeV

const double FWeak::gIn_en=6.671666e-2;///  20.015% / 3
const double FWeak::gIn_mn=6.671666e-2;///  1/15
const double FWeak::gIn_tn=6.671666e-2;///  1/15

const double FWeak::gOut_e=0.03367;
const double FWeak::gOut_m=0.03367;
const double FWeak::gOut_t=0.03371;
const double FWeak::gOut_q=0.6989;//sum Gamma_out_quark_i (hadron part)

//class FWeak - secondary quontities and methods

double FWeak::G2=0;//G^2
double FWeak::Mz2=0;//Mz^2
double FWeak::Mw2=0;//Mw^2
double FWeak::GammaZ2=0;//GammaZ^2
double FWeak::Cz0=0;//GammaZ^2*Mz^2
double FWeak::Cz1=0;//GammaZ^2*G^2/(4Pi)*Mz^2
double FWeak::Cz4=0;//6Pi/Mz^2*GammaZ^2
double FWeak::Cz5=0;//GammaZ^2/Mz^2
double FWeak::Cw0=0;// 1/Mw^2
double FWeak::Cw1=0;//G^2/(4Pi)*Mw^2
double FWeak::Cw2=0;//2.0*G2/Pi*Mw2*Mw2
double FWeak::Cw3=0;//2.0*G^2/Pi/Mw^2


int FWeak::isInit=0;
void FWeak::Init()
{
	if(isInit)
		return;
	isInit=1;
	G*=(units.Eunit*units.Eunit);
	Mz/=units.Eunit;
	Mw/=units.Eunit;
	GammaZ/=units.Eunit;
	G2=G*G;
	
	Mz2=Mz*Mz;
	Mw2=Mw*Mw;
	GammaZ2=GammaZ*GammaZ;
	Cz0=GammaZ2*Mz2;
	Cz1=GammaZ2*G2/4.0/Pi*Mz2;
	Cz4=6.*Pi/Mz2*GammaZ2;
	Cz5=GammaZ2/Mz2;
	Cw0=1.0/Mw2;
	Cw1=G2/4.0/Pi*Mw2*Mw2;
	Cw2=2.0*G2/Pi*Mw2*Mw2;
	Cw3=2.0*G2/Pi/Mw2;
}
//////// end FWeak
/***************************************************************************/

const double EMConstants::alpha = 0.00729735307639648;
// C1=0.5*Pi*alpha^2=3/16*TompsCS*Me^2
const double  EMConstants::C1=8.3647043703e-5;

EMConstants::EMConstants():
Me(ParticleData::getParticleMass(EElectron)),
Me2(Me*Me),
Mmu(105.658389/units.Eunit),
Mmu2(Mmu*Mmu)
{

}
