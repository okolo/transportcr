#include "PhotoDisintegration.h"
#include "Nucleus.h"
#include "Ranges.h"
#include "Medium.h"

namespace couplings{

PhotoDisintegration::PhotoDisintegration(void) :
fPDMap(*(new CPhotoDisintegrationMap()))
{
	const CParticleList& pf = *CParticleList::Instance();
	for(int secondary = ENeutron; secondary<EEndNuclei-1; secondary++)
	{
		if(!pf.IsEnabled((TParticle)secondary))
			continue;
                int no = 2 + secondary - EStartNuclei;
		const CPDMSecLine* canals = fPDMap.getIncome(no);
		int noOfCanals = canals ? canals->length() : 0;
		for(int i=0; i<noOfCanals; i++)
		{
			const TPhotoDisintegrationMapEntry& entry = (*canals)(i);
			const CPhotoDisintegrationCanal& canal = *(entry.iCanal);
			TParticle primary = canal.getParticle();
			if(pf.IsEnabled(primary))
				AddChannel(new Channel_A_x(this, primary, (TParticle)secondary, entry));
		}
	}
	for(int primary = EStartNuclei; primary<EEndNuclei; primary++)
	{
		if(!pf.IsEnabled((TParticle)primary))
			continue;
		AddChannel(new Channel_A_A(this, (TParticle)primary));
	}
}

void PhotoDisintegration::Channel_A_x::Coef(CMatrixAddOnlyView& aCoef) const
{
	const CPhotoDisintegrationCanal& canal = *(fEntry.iCanal);
	double rate = fEntry.iRate;
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		aCoef.Add(iPrim, iPrim, canal.R(iPrim)*rate);
	}
}

void PhotoDisintegration::Channel_A_A::Coef(CMatrixAddOnlyView& aCoef) const
{
	const CPDMLine* canals = fCoupling.fPDMap.getOutcome(CNucleus::getA(Primary()));
	int noOfCanals = canals?canals->length():0;
	if(!noOfCanals)
		return;
			
	const int nn = Ranges().nE();
	for(int iPrim = Ranges().nMinInteractionN(); iPrim<nn; iPrim++)
	{
		double result = 0.;
		for(int i=0; i<noOfCanals; i++)
			result += (*canals)(i).R(iPrim);
		aCoef.Add(iPrim, iPrim, -result);
	}
}

PhotoDisintegration::~PhotoDisintegration(void)
{
	delete &fPDMap;
}

void PhotoDisintegration::SetBackgrounds(const Medium& aPropagCoef)
{
	Coupling::SetBackgrounds(aPropagCoef);
	fPDMap.recalculate(&aPropagCoef.BackgroundIntegral());
}

}//namespace couplings{
