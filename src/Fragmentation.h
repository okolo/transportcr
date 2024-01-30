#if !defined(FRAGMENTATION_H_INCLUDED)
#define FRAGMENTATION_H_INCLUDED

class CFragmentation
{
public:
	CFragmentation(double _Ejet/*MeV*/,double _a,double _b);
	virtual double F(double x);

////atributes
protected:
	CFragmentation(){;};
	virtual void Normalize(double _c=1.0);
	static const double Nc;//number of quark colors
	static const double Nf;//number of quark flavors
	static const double Leff_MeV;//L_eff in MeV
	double norm,ksiM,sigma,sigma2;

};

class CFragmentationQCD : public CFragmentation
{
public:
    CFragmentationQCD(double _Ejet/*MeV*/);
};


class CFragmentationSUSY : public CFragmentation
{
public:
    CFragmentationSUSY(double _Ejet/*MeV*/);
};

class CFragmentationHill : public CFragmentation
{
public:
    CFragmentationHill(){;};
	virtual double F(double _x);
};
#endif //#if !defined(FRAGMENTATION_H_INCLUDED)
