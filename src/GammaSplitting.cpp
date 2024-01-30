/* 
 * File:   GammaSplitting.cpp
 * Author: ok
 * 
 * Created on October 24, 2011, 12:17 PM
 */

#include "GammaSplitting.h"
#include "Units.h"
#include "Parameters.h"

namespace couplings
{
    
double GammaSplitting::E0_eV = 1e19;
double GammaSplitting::L0_Mpc = 30;

GammaSplitting::GammaSplitting() {
	Reader()->readDoublePar("GammaSplittingE0", E0_eV);
	Reader()->readDoublePar("GammaSplittingL0", L0_Mpc);
    AddChannel(new Channel_gamma_gamma(this));
}

GammaSplitting::~GammaSplitting() {
}

void GammaSplitting::Channel_gamma_gamma::Coef(CMatrixAddOnlyView& aCoef) const
{
	const int nn = Ranges().nE();
        double E0 = E0_eV*1e-6/units.Eunit;
        double L0 = L0_Mpc*units.Mpc_cm /units.Lunit;
        const CBinning& E = Ranges().midE();
        //double binDifD = BC().s*log10(3.);
        int binDif = (int)(BC().s*log10(3.)+0.5);
        

	for(int iPrim = Ranges().nMinInteraction(); iPrim<nn; iPrim++)
	{
            double decayRate = E[iPrim]/E0;
            decayRate = decayRate*decayRate*decayRate/L0;
            aCoef.Add(iPrim, iPrim, -decayRate);
            if(iPrim>=binDif)
                aCoef.Add(iPrim-binDif, iPrim, 3.*decayRate);
	}
}

}//end namespace couplings


