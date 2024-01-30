// Stecker98IROSpectrum.cpp: implementation of the CStecker98IROSpectrum class.
//
//////////////////////////////////////////////////////////////////////

#include "Stecker98IROSpectrum.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#define NUMBER_OF_POINTS 63

static const double xArray[] = {
9.574852126495769E11	,
1.0840278284471202E12	,
1.2146647382177715E12	,
1.3610447882903398E12	,
1.5731332071794941E12	,
1.8371767615786367E12	,
2.058576146849856E12	,
2.3793592422370005E12	,
2.5317113923103096E12	,
2.721828371698936E12	,
3.1786754043900566E12	,
3.4888110431171777E12	,
3.829205862846499E12	,
4.335279189003681E12	,
4.758261849795373E12	,
5.387120544601277E12	,
6.35674077335715E12		,
7.34729678702299E12		,
8.064154689158595E12	,
9.12992491211411E12		,
1.0124902169134488E13	,
1.3112935051592441E13	,
1.8258097654881934E13	,
2.2454468786842566E13	,
3.030967724478728E13	,
3.766352243237016E13	,
4.3532540844697086E13	,
5.031611463845671E13	,
5.755828121517693E13	,
6.252399512835922E13	,
6.933784669555311E13	,
7.610297198742447E13	,
7.93178747531347E13		,
9.263103442812088E13	,
1.163023525170387E14	,
1.3031802585417752E14	,
1.4452005633666872E14	,
1.6026982105242194E14	,
1.68777118092752E14		,
1.852442915732296E14	,
2.0972644221905744E14	,
2.4747480875308475E14	,
2.660587450018289E14	,
3.172106127809051E14	,
3.5913362127203256E14	,
4.10824923894294E14		,
4.797800973314469E14	,
5.661350403618632E14	,
6.543545366055695E14	,
7.256659535359305E14	,
8.562775450238892E14	,
9.594681176756681E14	,
1.0975675711646411E15	,
1.2298361642752315E15	,
1.3638634454515985E15	,
1.528223530301249E15	,
1.6947690025729415E15	,
1.8794645646608412E15	,
2.1059600810189875E15	,
2.264105583121814E15	,
2.5108473235396365E15	,
2.671618753703095E15	,
2.872241972082083E15	};

static const double yArray[] = {
2.274694496808945E-9,
2.793533257985222E-9,
3.358565682996543E-9,
4.009377262024964E-9,
4.994149495797963E-9,
5.36081437182297E-9,
5.633382376146457E-9,
5.554121610298015E-9,
5.398929946034567E-9,
5.101434415464322E-9,
4.458931178088874E-9,
3.981071705534974E-9,
3.579696061042172E-9,
3.10675990983383E-9,
2.7154771573268926E-9,
2.2426899184136027E-9,
1.8004647182119289E-9,
1.5406085639499548E-9,
1.4251026703029986E-9,
1.3465756882423677E-9,
1.2997090825113858E-9,
1.3951322711861721E-9,
1.5515623437261944E-9,
1.6420433388871625E-9,
1.7501566522976255E-9,
1.8522188841040355E-9,
1.960232982646278E-9,
2.1191116846058047E-9,
2.3734747473159386E-9,
2.6209669129974093E-9,
3.084826686197942E-9,
3.579696061042172E-9,
3.656595522033888E-9,
4.458931178088874E-9,
6.004288319960372E-9,
6.9183097091893634E-9,
7.585775750291836E-9,
8.142714930301282E-9,
8.376776400682916E-9,
8.67883713269597E-9,
9.055721899982202E-9,
9.448973160341253E-9,
9.382264934578687E-9,
9.055721899982202E-9,
8.74054396263963E-9,
8.085228676856247E-9,
7.585775750291836E-9,
7.270067807853826E-9,
6.9183097091893634E-9,
6.820970226563626E-9,
7.01703828670383E-9,
7.117175782354046E-9,
6.869467559027169E-9,
6.309573444801936E-9,
5.836518227320602E-9,
4.95889160624696E-9,
4.124626382901353E-9,
3.2187875118212363E-9,
2.459060596250654E-9,
1.9326528261662556E-9,
1.3465756882423677E-9,
9.448973160341251E-10,
5.713774242363902E-10	};


CStecker98IROSpectrum::CStecker98IROSpectrum(double aZmax):
m_Zmax(aZmax),
m_X(NUMBER_OF_POINTS),
m_Y(NUMBER_OF_POINTS),
m_spectrum(NULL)
{
}

bool CStecker98IROSpectrum::init()
{
	m_X.copy(xArray);
	m_Y.copy(yArray);

	m_X *= (1./(241803*1e9));//Hz to eV

	for(int i=0; i<m_X.length(); i++)
		m_Y[i] *= (256153.72/m_X[i]);
	/// converting nu * f(nu) to n(k) * k
	/// [nu * f(nu)] = W/m^2/sr
	/// the coeffitient 256153.72 is equal to 4Pi/(speedOfLight/(m/sec))/1.6e-19/(m^3/cm^3)
	/// todo: the number 256153.72 is not very accurate and should be replaced by exact value

	m_spectrum = new CDefaultTableFunc(m_X, m_Y);
	return true;
}

CStecker98IROSpectrum::~CStecker98IROSpectrum()
{
	delete m_spectrum;
}

double CStecker98IROSpectrum::F(double E, double z)
{
	if(z>m_Zmax)
		return 0.;
	return m_spectrum->f(E);
}

//returns maximal background energy in eV
double CStecker98IROSpectrum::MaxE(double aZmax) const
{
	return m_spectrum->lastX();
}

//returns minimal background energy in eV
double CStecker98IROSpectrum::MinE(double aZmax) const
{
	return m_spectrum->firstX();
}
