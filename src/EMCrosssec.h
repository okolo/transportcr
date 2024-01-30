#ifndef EM_CROSSEC_H_INCLUDED
#define EM_CROSSEC_H_INCLUDED

class EMCrosssec
{
public:
	EMCrosssec(double aMass);
	double RPP(double E,double b);
	double PPPdrIntegral(double E, double b, double r1, double r2);
	double lowerPPPdrIntegral(double E, double b, double r);
	double lowerPPPrdrIntegral(double E, double b, double r);
	double PPPe(double E,double b,double nE);
	double rMaxICSph(double Eb);
	double PICSph(double E,double b,double nE);
	double EnergyLossICS(double Eb);
	double PICSe(double E,double b,double nE);
	double RICS(double Eb);
	double RICS(double E,double b);
	double intPrdrICSph(double E, double b);
	double PICSphR(double Eb, double r);
	double PICSeR(double Eb, double r);
private:
	const double fMass;
	const double fSquareMass;
    /* C1=0.5*Pi/137.035989561/137.035989561=
       =3/16*TompsCS*Me^2=Pi/2*alpha^2 */
	static const double  C1;
	static const double ln2;
	double li2(double aX);
	double intPICSph_Jones(double A, double q);
	double intPICSph_Jones(double E, double b, double nE);//return intergal of P from 0 to nE
	//maximal nE/E for secondary photon in ICS (see PICSphR(double Eb, double r))
	double intPrdrICSph_Jones(double E, double b);
	double PPPe_old(double E,double b,double nE);
	double PPPeAharonyan(double E,double b,double nE);
	double RPP_Gould(double E,double b);
	double intPPPdr(double r, double W);
	double intPPPrdr(double r, double W);
	double RICS_JonesSeriesZero(double A);
	double RICS_JonesSeriesInf(double A);
	double F_ics(double b, double _m);
	double F_I_PICSe(double x,double k);
	double I_PICSe(double x_min,double x_max,const double s);
	double I_PICSph(double x_min,double x_max,const double s);
	double Pe(double _x,double _s);
	double Pph(double _x,double _s);
	double Re(double _s);
	//double F_ics(double b, double m/*needed if 1-b (if it is small)*/);
	double RPP_old(double E,double b);
};

#endif //end of file
