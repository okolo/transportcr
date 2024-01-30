#if !defined(MUPP_H_INCLUDED)
#define MUPP_H_INCLUDED

//#include "Prodspec.h"
#include "const.h"

class CMuPairProduction : EMConstants
{
enum TProducts{
	en=0,
	mn=1,
	e=1//same spectrum
	};
protected:
	double** coef[2];
	int s_min;
	int s_max;
	double GetP(int _nE,int _E,int _b,TProducts product);
	double PmuPP(double x,double s);
public:
	CMuPairProduction(const double aMinK, const int nK);
	~CMuPairProduction();
	inline double GetPen(int _nE,int _E,int _b){return GetP(_nE,_E,_b,en);};
	inline double GetPmn(int _nE,int _E,int _b){return GetP(_nE,_E,_b,mn);};
	inline double GetPe(int _nE,int _E,int _b){return GetP(_nE,_E,_b,e);};
};

#endif //end of file
