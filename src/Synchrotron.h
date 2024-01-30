#pragma once
#include "Coupling.h"
#include "Addfunc.h"
#include "const.h"
#include "TableFunc.h"

namespace couplings
{

class Synchrotron : public Coupling, protected CouplingParameters, EMConstants
{
	class Channel_e_e : public CouplingChannelT<Synchrotron>
	{
	public:
		Channel_e_e(Synchrotron* aCoupling, TParticle aPrim) : CouplingChannelT<Synchrotron>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_e_gamma : public CouplingChannelT<Synchrotron>
	{
	public:
		Channel_e_gamma(Synchrotron* aCoupling, TParticle aPrim) : CouplingChannelT<Synchrotron>(aCoupling, aPrim, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	Synchrotron();
    inline double SynchroEnergy(double Ee, double B)
    /* Synchrotron radiation characteristic energy
     * Ee - electron energy, B - random magnetic field amplitude */
    {
        const double c1=0.128136819178567;//3/2*sqrt(alpha)
        return c1*B/Me*Ee/Me2*Ee;
    }
	inline double SynchroSpectrum(double Ee, double Eph, double B)
    {
        /* Synchrotron radiation spectrum dP(Ee,Eph)/dEph
         * Ee - electron energy, Eph - photon energy, B - random magnetic field amplitude */
        const double c2=0.000171841794354763;//3^(1/2)/2/Pi*alpha^(3/2)
        double Ec=SynchroEnergy(Ee, B);
        return c2*B/Me*fG(Eph/Ec);
    }


	void SetBackgrounds(const Medium& aPropagCoef);
private:
//	class Data
//	{
//		int ac;
//		double kx;
//		double k;
//		double *x;
//		double *y;
//		double ymin;
//		double ymax;
//		double lbackground;
//		double rbackground;
//	public:
//		Data():ac(0),kx(0),k(0),x(0),y(0),ymin(0),ymax(0),lbackground(0),rbackground(0){};
//		inline bool IsInitialized() const { return x!=0;}
//		int InitData(const char *_fileName,double lbackground=0,double rbackground=0);
//		int InitDataL(char *_fileName,double _lbackground);
//		int InitDataR(char *_fileName,double _rbackground);
//		int InitDataC(char *_fileName);
//		virtual double GetData(double x);
//		virtual ~Data();
//	};

//  static double G(double x);
	static const double B_const;
	double 	fPrevB;
//	static Data fGData;
    CLinearFunc fG;
    SafePtr<CLinearFunc> fIntG;

	CMatrix Syn;//synchrotron radiation coef for photons
	CMatrix Syn_e;//synchrotron radiation coef for electrons
	CVector aeSyn;//synchrotron radiation coef for electrons
};
}
