/*
 * Inoue12IROSpectrum.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: ok
 */

#include "Inoue12IROSpectrum.h"

Inoue12IROSpectrum::Inoue12IROSpectrum(std::string aDataFile):
MatrixBackground(DATA_DIR + aDataFile, false, true)
{
}

double Inoue12IROSpectrum::F(double E, double z)
{
	return MatrixBackground::F(E,z)/E;
}

Inoue12BaselineIROSpectrum::Inoue12BaselineIROSpectrum():
Inoue12IROSpectrum("EBL_Inoue/EBL_proper_baseline.dat")
{
}

Inoue12LowPop3IROSpectrum::Inoue12LowPop3IROSpectrum():
Inoue12IROSpectrum("EBL_Inoue/EBL_proper_low_pop3.dat")
{
}

Inoue12UpperPop3IROSpectrum::Inoue12UpperPop3IROSpectrum():
Inoue12IROSpectrum("EBL_Inoue/EBL_proper_up_pop3.dat")
{
}


