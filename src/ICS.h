#pragma once
#include "Coupling.h"
#include "TableFunc.h"
#include "EMCrosssec.h"
#include "const.h"

class CBackgroundTable;

namespace couplings
{
class ICSold : public Coupling, protected CouplingParameters, EMConstants
{
	class Channel_e_e : public CouplingChannelT<ICSold>
	{
	public:
		Channel_e_e(ICSold* aCoupling, TParticle aPrim) : CouplingChannelT<ICSold>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_e_gamma : public CouplingChannelT<ICSold>
	{
	public:
		Channel_e_gamma(ICSold* aCoupling, TParticle aPrim) : CouplingChannelT<ICSold>(aCoupling, aPrim, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	ICSold(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	virtual void DebugOutput(const char* aFolder) const;
private:
	void InitCoef(const CBackgroundTable& aBackgr);
	void AdjustCoef();//ICS part of former PropagCoef::LeeAdjustment()
	void FinalCoefAdjustment();//part from PropagCoef::coefficients()
	CVector	f_a;
	CVector	f_intPrdrGamma;
	CMatrix	f_be;
	CMatrix	f_bph;
	EMCrosssec cs;
};

class ICScel : public Coupling, protected CouplingParameters, EMConstants
{
	class Channel_e_e : public CouplingChannelT<ICScel>
	{
	public:
		Channel_e_e(ICScel* aCoupling, TParticle aPrim) : CouplingChannelT<ICScel>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_e_gamma : public CouplingChannelT<ICScel>
	{
	public:
		Channel_e_gamma(ICScel* aCoupling, TParticle aPrim) : CouplingChannelT<ICScel>(aCoupling, aPrim, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class PFunction : public IFunction
	{
	public:
		double (* P) (double E, double b,double nE);
		double E;
		double b;
	protected:
		virtual double f(double _x) const
		{
			return P(E, b, _x);
		}
	};
public:
	ICScel(void);
	virtual void SetBackgrounds(const Medium& aPropagCoef);
	virtual void DebugOutput(const char* aFolder) const;
private:
	void InitCoef(const CBackgroundTable& aBackgr);
	void AdjustCoef();//ICS part of former PropagCoef::LeeAdjustment()
	static double PICSe(double E,double b,double nE);
	static double PICSph(double E,double b,double nE);

	CVector fEnergyLossRate;//-1/E dE/dt for electron
	CVector	f_a;
	CVector	f_intPrdrGamma;
	CMatrix	f_be;
	CMatrix	f_bph;
	EMCrosssec cs;
	static double ICS_CEL_MAX_BINS_DIFF;
};

class ICSsplitted : public Coupling, protected CouplingParameters, EMConstants
{
	class Channel_e_e : public CouplingChannelT<ICSsplitted>
	{
	public:
		Channel_e_e(ICSsplitted* aCoupling, TParticle aPrim) : CouplingChannelT<ICSsplitted>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class Channel_e_gamma : public CouplingChannelT<ICSsplitted>
	{
	public:
		Channel_e_gamma(ICSsplitted* aCoupling, TParticle aPrim) : CouplingChannelT<ICSsplitted>(aCoupling, aPrim, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};

	class AKern : public IFunction
	{
	public:
		AKern(IFunction* aR, double aE):fR(aR),fE(aE){}
		double f(double b) const{return (*fR)(fE*b)/(b*BC().sF);}
	private:
		IFunction*	fR;
		double		fE;
	};

	class RhoKern : public IFunction
	{
	public:
		RhoKern(IFunction* aR, IFunction* aRho, double aE):fR(aR),fRho(aRho),fE(aE){}
		double f(double b) const{return ((*fR)(fE*b))*((*fRho)(fE*b))/(b*BC().sF);}
	private:
		IFunction*	fR;
		IFunction*	fRho;
		double		fE;
	};

	class RKern : public IFunction
	{
	public:
		RKern(double aP(double,double), double aW):fP(aP),fW(aW){}
		double f(double r) const{return fP(fW, r);}
		inline void SetW(double aW) { fW = aW; }
	private:
		double (*fP)(double,double);
		double fW;
	};

public:
	ICSsplitted();
	~ICSsplitted();
	virtual void SetBackgrounds(const Medium& aPropagCoef);
private:
	void InitCoef(const CBackgroundTable& aBackgr);
	void Init(double aKmin, double aKmax);
	static double PICSphR(double Eb, double r);
	static double PICSeR(double Eb, double r);
	CVector		 fEnergyLossRate;//-1/E dE/dt for electron
	CLinearFunc* fRcel;
	CLinearFunc* fRdiscr;
	CLinearFunc* fRho;
	CVector		 fRcelV;
	CVector		 fRdiscrV;
	CVector		 fW;
	CVector		 fRhoV;
	CVector		 fa_discr;
	CMatrix		 fb_discr;
	CMatrix		 f_bph;
	EMCrosssec	 cs;
	static double ICS_CEL_SPLITTED_MAX_BINS_DIFF;
};

}//namespace couplings
