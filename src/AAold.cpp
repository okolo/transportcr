//
// Created by ok on 14/10/2015.
//

#include "AAold.h"
#include "Units.h"
#include "Nucleus.h"
#include "Ranges.h"
#include "Medium.h"

#ifdef AA_FORTRAN

#define Sighad_fortran sighad_
#define Sigtap2_fortran sigtap2_

extern "C" void Sigtap2_fortran(int*);
extern "C" double Sighad_fortran(int*,double*,double*,double*,double*,double*); // defined in crn6.f

#else

void Sigtap2_fortran(int*)
{
    throw "Fatal error in Sigtap2_fortran: Fortran integration is disabled";
}

double Sighad_fortran(int*,double*,double*,double*,double*,double*)
{
    throw "Fatal error in Sighad_fortran: Fortran integration is disabled";
    return 0.;
}

#endif

// Barashenkov & Polanski pA total cross section  IMOS20020502
double sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E)
{
    if(TA<4)
    {//sighad_cc can only calculate cross section for nuclei with A>=4
        TA=4;
        TZ=2;
    }
    return( Sighad_fortran(&IS, &PA, &PZ, &TA, &TZ, &E) );
}

// initialization of the Barashenkov & Polanski cross section code
void sigtap_cc(int ISS)
{
    Sigtap2_fortran(&ISS);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// from nucleon_cs.cc *                               galprop package * 2001/05/11
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// parametrization of the pp-, pA-, AA-, (anti_p)p-, and (anti_p)A total
// inelastic cross sections, as well as (anti_p) annihilation cross sect.
// A sample driver program is in the bottom of this file.
//    INPUT:
// option -pA total inelastic cross section
//       =0 -[L83];
//       =1 -[WA96] for Zt>5 and [BP01] for Zt<=5;
//       =2 -[BP01];
// Ek -kinetic energy per nucleon of the beam momentum particles (GeV);
// Zp=+/-1 is the charge of beam and At is the atomic number of target nuclei
// (Zp=At=1 for pp-collisions, Zp= -1 for antiprotons);
//    OUTPUT:
// PP_inel  [mbarn] is the proton-proton inelastic cross sect.;
// PA_inel  [mbarn] is the proton-nucleus inelastic cross sect. (p+A2);
// aPP_non  [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.);
// aPA_non  [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.);
// aPP_ann  [mbarn] is the antiproton-proton annihilation cross sect.,
// aPA_ann  [mbarn] is the antiproton-nucleus annihilation cross sect.,
//    REFERENCES:
// [W79] Westfall et al. 1979, PRC, 19, 1309
// [L83] Letaw et al. 1983, ApJS, 51, 271
// [TN83] Tan & Ng 1983, J.Phys.G:Nucl.Phys., 9, 227
// [PDG00] D.E.Groom et al. 2000, Europ. Phys. J. C15, 1
// [MO97] Moiseev & Ormes 1997, Astropart. Phys.6, 379
// [WA96] Wellisch H.P., Axen D. 1996, PRC 54, 1329; Wellish 1999, private comm.
//    (typo corrections & code)
// [BP01] V.S.Barashenkov, A.Polanski code used for calc. of the pA total inelastic
//    cross sections
void nucleon_cs(int option, double Ek, int Zp, int Zt, int At,
                double *PP_inel,double *PA_inel,double *aPP_non,double *aPA_non,double *aPP_ann,double *aPA_ann)
{
    const double Mp   = 0.938;               // GeV, Proton rest mass
    *PP_inel= *PA_inel= *aPP_non= *aPA_non= *aPP_ann= *aPA_ann=0.;
    if(Ek <= 0.) return;

    double U,C1,s2,Z,A,b0,Fcorr,rN,s0,p1,p2,p3,p4,p5,x,f1,f2;
    double aPP_inel=0.;
    //double aPA_inel=0.;
    //double Cp;
    double Emev=Ek;
    double PZ=fabs(1.*Zp),Em=1000.*Ek, TZ=Zt, TA=At;
    int ISS=2;

//## Proton-Proton INELASTIC cross section, mb [TN83,p.234]
    if(Ek > 0.3)
    {
        U = log((Ek+Mp)/200.);
        *PP_inel = 32.2*(1.+0.0273*U);
        if(U >= 0.) *PP_inel += 32.2*0.01*U*U;
        if(Ek < 3.) *PP_inel /= 1.+0.00262/pow(Ek,17.9+13.8*log(Ek)+4.41*pow(log(Ek),2));
    }
    if(Zp*At == 1) return;

//## Proton-Nucleus INELASTIC cross section, mb
    switch(option)
    {
        case 0:                                        // [L83]
            C1 = (At == 4) ? 0.8 : 1.;                  // a correction for He
            if(At == 9) C1 = 1.+0.75*exp(-Ek/0.075);    // a correction for Be
            *PA_inel = C1 *45. *pow(TA,0.7) *(1.+0.016*sin(5.3-2.63*log(TA)));
            if(Ek < 5.) *PA_inel *= 1.-0.62 *exp(-Ek/0.2) *sin(10.9/pow(Ek*1.e3,0.28));
            if(At == 4) *PA_inel = (Ek > 0.01) ?        // pHe, my fit
                                   111.*(1.-(1.-sin(9.72*pow(log10(Ek*1000.),0.319)-4.14))*exp(-3.84*(Ek-0.1))) : 0.;
            break;

        case 1:                                        // Zt>5 [WA96], Zt<=5 [BP01]
            if(Zt>5)
            {
                b0 = 2.247-0.915*(1.-pow(TA,-1./3.));
                Fcorr = (1.-0.15*exp(-Emev))/(1.-0.0007*At);   // high-energy correction
                rN = (At-Zt>1.5) ? log(TA-Zt) : 1.;
                s0 = Pi*10.*pow(1.36,2.)*Fcorr*rN*(1.+pow(TA,1./3.)-b0*(1.-pow(TA,-1./3.)));
                p1 = 8.-8./At-0.008*At;
                p2 = 2.*(1.17-2.7/At-0.0014*At);
                p3 = 0.8+18./At-0.002*At;
                p4 = 5.6-0.016*At;
                p5 = 1.37*(1.+1./At);
                x = log10(Emev);
                f1 = 1./(1.+exp(-p1*(x+p2)));                 // low-energy return to zero
                f2 = 1. +p3 *( 1. -1./(1.+exp(-p4*(x+p5))) ); // low-energy rise
                *PA_inel = f1*f2*s0;
                break;
            }

        case 2:                                        // [BP01]
        default:
            if (Em<14.)  Em=14.;
            if (Em>1.e6) Em=1.e6;
            *PA_inel = sighad_cc(ISS,PZ,PZ,TA,TZ,Em);   // IMOS20020502
    }
    if(Zp*At >= 1) return;

//## AntiProton-Proton ANNIHILATION cross section [TN83]
    if(Ek < 10.) *aPP_ann = 661.*(1.+0.0115/pow(Ek,0.774)-0.948*pow(Ek,0.0151)); // 0.1GeV<Ek<12GeV
    else
    {
// assuming aPP_ann = aPP_tot -PP_tot (i.e., aPP_elast = PP_elast); (aPP_tot-PP_tot) from [PDG00]
        s2 = 2.*Mp*(Ek+2*Mp);                   // square of the total CMS energy
        *aPP_ann = 2*35.43/pow(s2,0.560);
    }

//## AntiProton-Proton TOTAL INELASTIC cross section
    aPP_inel = *PP_inel + *aPP_ann;
    if(Ek <= 14.)
    {
        aPP_inel = 24.7*(1.+0.584/pow(Ek,0.115)+0.856/pow(Ek,0.566));
        if(*aPP_ann > aPP_inel) *aPP_ann = aPP_inel;
    }

//## AntiProton-Proton TOTAL INELASTIC NON-ANNIHILATION cross section
    *aPP_non = aPP_inel - *aPP_ann;
    if(*aPP_non < 0.) *aPP_non = 0.;

//## AntiProton-NUCLEUS cross sections
    if(At > 1)
    {
//## AntiProton-NUCLEUS TOTAL INELASTIC NON-ANNIHILATION cross section
        *aPA_non = *PA_inel;

//## AntiProton-NUCLEUS ANNIHILATION cross section on 12C-nucleus [mb] (0.4<Pp<300) [MO97]
        A = At;
        Z = Zt;                               // Z = 0.59*pow(A,.927);  for Z > 2 nuclei
        if(At == 4) { Z = 2.; A = 3.30; }     // Modified to agree with HE p-He cs / imos
        *aPA_ann = pow(A,2./3.)               // Scaling to other nuclei
                   //         *(48.2 +19./pow(Ek-0.02,0.55)
                   *(48.2 +19./pow(Ek,0.55)           // modified to agree w. He@<100 MeV / imos
                     +(0.1-0.18/pow(Ek,1.2))*Z +0.0012/pow(Ek,1.5)*Z*Z)  - *aPA_non;
        if(*aPA_ann < 0.) *aPA_ann = 0.;
        if(*aPA_ann < *aPP_ann) *aPA_ann = *aPP_ann;
        if(At == 4 && Ek > 5.)  *aPA_ann = *aPP_ann;
/*
//## my fit to AntiProton-NUCLEUS total cross section on 12C-nucleus [mb] (0.4<Pp<300)
         double Pp =sqrt(pow(Ek+Mp,2)-Mp*Mp); // GeV,kin. momentum per nucleon
         *aPA_ann = (Pp > 40.) ? 236.*(1.+6.9e-2*exp(-Pp/100.)) :
            209.7*pow(.542/Pp,.565)+29.6*cos(log10(Pp/9.29)*5.11)+257.9;
         *aPA_ann *= pow(At/12.,2./3.);                             // scaling to other nuclei
*/
    }
    return;
}



// Ep - energy per nucleon E/A
// (sighad_cc is actually calculated in the rest frame of nucleus
// and so Ep = E/A should be passed to sighad_cc)
//
// A - nucleus atomic number
// Z - nucleus charge
double A_p_cs(double Ep, int A, int Z)
{
    Ep*=units.Eunit;//converting to MeV
//PROJECTILE PARTICLE KINETIC ENERGY in sighad (14MEV < Ep < 1TEV)
    if(Ep<14)
        Ep=14;//using left value
    else if(Ep>1e6)
        Ep=1e6;//using right value

    double result = sighad_cc(
            2,//inelastisic cross section only
            1,//proton  (A1=1.0 FOR NUCLEON OR 0<A1<0.2 FOR PIONS)
            1,//proton  (Z1-PROJECTILE CHARGE NUMBER)
            A,
            Z,
            Ep
    );

    if(result>0)
        return result;
    return 0.;
}

double A_p_csVerGalprop(double Ep, int A, int Z)
{
    int option = 2;//V.S.Barashenkov, A.Polanski code used for calc. of the pA total inelastic cross sections
    double Ek = Ep*units.Eunit*1e-3;//converting to GeV
    int Zp=1;//proton
    double PP_inel, PA_inel, aPP_non, aPA_non, aPP_ann, aPA_ann;

    nucleon_cs(option, Ek, Zp, Z, A, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);

    if(PA_inel>0)
        return PA_inel;//breakpoint
    return 0.;
}

namespace couplings {

    bool AAold::fInitialized = false;

    AAold::AAold() {
        Init();
        const CParticleList &pf = *CParticleList::Instance();
        for (int primary = EStartNuclei; primary < EEndNuclei; primary++) {
            if (!pf.IsEnabled((TParticle) primary))
                continue;
            AddChannel(new Channel_A_A(this, (TParticle) primary));
        }
    }

    void AAold::Init() {
        if (fInitialized)
            return;
        //ensure that data file exists
        FILE *f = FopenL("tables/barpol.dat", "rt");//don't use DATA_DIR macro here since it is not defined in crn6.f
        //FILE* f = FopenL("barpol.dat", "rt");//don't use DATA_DIR macro here since it is not defined in crn6.f
        fclose(f);
        sigtap_cc(-1);
        fInitialized = true;
    }

    void AAold::Channel_A_A::Coef(CMatrixAddOnlyView &aCoef) const {
        int A = CNucleus::getA(Primary());
        int Z = CNucleus::getZ(Primary());
        double mult = (*Background()).protonsConc * units.barn *
                      1e-3;//converting mbarn to internal units and multiplying by concentration
        const int nn = Ranges().nE();
        for (int iPrim = Ranges().nMinInteractionN(); iPrim < nn; iPrim++) {
            double E = Ranges().midE()[iPrim];
            //double suppression = mult*A_p_cs(E, A, Z);
            double suppression = mult * A_p_csVerGalprop(E, A, Z);
            aCoef.Add(iPrim, iPrim, -suppression);
        }
    }
}