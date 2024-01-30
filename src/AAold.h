//
// Created by ok on 14/10/2015.
//

#ifndef PROPAGATION_AAOLD_H
#define PROPAGATION_AAOLD_H


#include "Coupling.h"

#ifdef NO_FORTRAN
#undef AA_FORTRAN
#endif
namespace couplings {
	class AAold : public Coupling {
		class Channel_A_A : public CouplingChannelT<AAold> {
		public:
			Channel_A_A(AAold *aCoupling, TParticle aPrim) : CouplingChannelT<AAold>(aCoupling, aPrim, aPrim) { };

			virtual void Coef(CMatrixAddOnlyView &aCoef) const;
		};

	public:
		AAold();

	private:
		static void Init();

		static bool fInitialized;
	};
}
#endif //PROPAGATION_AAOLD_H
