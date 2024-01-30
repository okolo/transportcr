/* 
 * File:   GammaSplitting.h
 * Author: ok
 *
 * Created on October 24, 2011, 12:17 PM
 */

#ifndef GAMMASPLITTING_H
#define	GAMMASPLITTING_H

#include "Coupling.h"
#include "TableFunc.h"

class CBackgroundTable;

namespace couplings
{

class GammaSplitting  : public Coupling, Parameters{
public:
    GammaSplitting();
    virtual ~GammaSplitting();
private:
    class Channel_gamma_gamma : public CouplingChannelT<GammaSplitting>
        {
	public:
                Channel_gamma_gamma(GammaSplitting* aCoupling) : CouplingChannelT<GammaSplitting>(aCoupling, EPhoton, EPhoton){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
    static double E0_eV;
    static double L0_Mpc;
};

}

#endif	/* GAMMASPLITTING_H */

