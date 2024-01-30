#pragma once
#include "Coupling.h"

class CPhotoDisintegrationMap;
struct TPhotoDisintegrationMapEntry;
namespace couplings{

class PhotoDisintegration : public Coupling
{
	class Channel_A_x : public CouplingChannelT<PhotoDisintegration>
	{
	public:
		Channel_A_x(PhotoDisintegration* aCoupling, TParticle aPrim, TParticle aSec,
			const TPhotoDisintegrationMapEntry& aEntry) : 
			CouplingChannelT<PhotoDisintegration>(aCoupling, aPrim, aSec),
				fEntry(aEntry){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	private:
		const TPhotoDisintegrationMapEntry& fEntry;
	};

	class Channel_A_A : public CouplingChannelT<PhotoDisintegration>
	{
	public:
		Channel_A_A(PhotoDisintegration* aCoupling, TParticle aPrim) : CouplingChannelT<PhotoDisintegration>(aCoupling, aPrim, aPrim){};
		virtual void Coef(CMatrixAddOnlyView& aCoef) const;
	};
public:
	PhotoDisintegration(void);
	~PhotoDisintegration();

	virtual void SetBackgrounds(const Medium& aPropagCoef);
private:
	CPhotoDisintegrationMap&	fPDMap;
};									 

}//namespace couplings{
