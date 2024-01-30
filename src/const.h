 /*
fundamental constants, included in F* classes
*/
#if !defined(CONST_H_INCLUDED_)
#define CONST_H_INCLUDED_ 

extern const double Pi;
class FWeak //weak interaction constants
{
	static int isInit;
public:
	static void InitGlobal();
	static double G;//Fermi constant / MeV^(-2)
	static const double gL_l;
	static const double gR_l;
	static const double gL_n;
	static const double gR_n;

	static double Mz;//MeV/
	static double Mw;//MeV
	static double GammaZ;//MeV
	
	// GammaZ_in/GammaZ and GammaZ_out/GammaZ for different particles
	// const names : gIn_#particle, gOut_#particle
	
	 
	static const double gIn_en;
	static const double gIn_mn;
	static const double gIn_tn;
	
	static const double gOut_e;
	static const double gOut_m;
	static const double gOut_t;
	static const double gOut_q;//sum Gamma_out_quark_i
	//secondary quontities
	
	static double G2;//G^2

	static double Mz2;//Mz^2
	static double Mw2;//Mw^2
	static double GammaZ2;//GammaZ^2
	
	static double Cz0;//GammaZ^2*Mz^2
	static double Cz1;//GammaZ^2*G^2/(4Pi)*Mz^2
	static double Cz4;//6Pi/Mz^2*GammaZ^2
	static double Cz5;//GammaZ^2/Mz^2
	
	static double Cw0;// 1/Mw^2
	static double Cw1;//G^2/(4Pi)*Mw^4
	static double Cw2;//2.0*G^2/Pi*Mw^4
	static double Cw3;//2.0*G^2/Pi/Mw^2
	
	////////////////////  methods
	
	static void Init();// converting to internal units system
	
};

class EMConstants
{
public:
	EMConstants();
	const double Me;
	const double Me2;
	const double Mmu;
	const double Mmu2;

	static const double alpha;
	// C1=0.5*Pi*alpha^2=3/16*TompsCS*Me^2
	static const double  C1;
};

#endif
