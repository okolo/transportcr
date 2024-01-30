// TimeZ.h: interface for the CTimeZ class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIMEZ_H__A1DF137A_F6C1_46E4_981A_5FDDDD0958B1__INCLUDED_)
#define AFX_TIMEZ_H__A1DF137A_F6C1_46E4_981A_5FDDDD0958B1__INCLUDED_

#include "debug.h"
#include "Addfunc.h"

class CLinearFunc;
class CVector;

class CTimeZ  
{
public:
	static void init(double _Lv, double _Hubble_in_km_s_Mpc);

// z -> t, dz -> dt and back transformation
	static double z2t(double _z);
	static double t2z(double _t);
	static double z2d(double _z);
	static double d2z(double _d);
	static double z2Dc(double _z);
	static double H(double _z);
	inline static double getLv(){return Lv;};
	inline static double getAgeOfUniverse(){return t0;};
	static void SuppressionTest();
	//returns energy loss rate -1/E dE/dt due to redshift at 
	inline static double eLossRate(double aZ);

	// -dt/dz
	static double diffT(double aZ);
	static double diffD(double a);
	//dDc/dz ( Dc is comoving distance as defined in astro-ph/9905116 )
	static double diffDc(double aZ);
protected:
	static void setLv(double Lv);
	static const int Accuracy;
	static double TestIntegralF(double r);
	static double TestIntegral(double r);
	static double Lv;// omega_lambda at z=0
	static double Lm;// omega_matter at z=0
	static double sqrtLv;//sqrt(Lv)
	static double sqrtLm;//sqrt(Lv)
	static double t0;//universe age at z=0
	static double Hubble;
	void initZ2D();
	void initZ2T();
	void initZ2Dc();
	SafePtr<CLinearFunc>	z2dFunc;
	SafePtr<CLinearFunc>	d2zFunc;
	SafePtr<CLinearFunc>	z2tFunc;
	SafePtr<CLinearFunc>	t2zFunc;
	SafePtr<CLinearFunc>	z2DcFunc;

	SafePtr<CVector>	zt;
	SafePtr<CVector>	tt;

	SafePtr<CVector>	z;
	SafePtr<CVector>	d;

	SafePtr<CVector>	zDc;
	SafePtr<CVector>	dDc;
};

inline double CTimeZ::eLossRate(double aZ)
{
	return 1./(1.+aZ)/diffT(aZ);
}

class Redshift
{
public:
	//uninitialized constructor
	Redshift();

	//constructor with redshift value z
	Redshift(double aZ);

	//redshift value
	inline double z() { ASSERT(m_z>=0.); return m_z;}

	//(1+z)^3
	inline double dv() { ASSERT(m_z>=0.); return m_dv;}

	//(1+z)
	inline double dl() { ASSERT(m_z>=0.); return m_dl;}

	//universe age at current redshift value z
	inline double t() { ASSERT(m_z>=0.); return m_t;}

	//set redshift value z
	void setZ(double aZ);

	//set redshift using value of universe age
	inline void setT(double aT) {setZ(CTimeZ::t2z(aT));}
private:
	double m_z;
	double m_dv;
	double m_dl;
	double m_t;
};

extern Redshift redshift;

#endif // !defined(AFX_TIMEZ_H__A1DF137A_F6C1_46E4_981A_5FDDDD0958B1__INCLUDED_)
