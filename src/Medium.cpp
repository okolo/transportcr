
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "Addfunc.h"
//#include "Prodspec.h"
#include <string.h>
#include "time.h"
#include "main.h"
#include "Concentrations.h"
#include "TimeZ.h"
#include "FilePtr.h"
#include "Background.h"
#include "Nucleus.h"
#include "ParticleList.h"
#include <fstream>
#include "SafeOutput.h"
#include "Log.h"
#include "PropagEngine.h"
#include "Parameters.h"
#include "Units.h"
#include "RadioBackground.h"
#include "Sarkar1005Background.h"
#include "PrimackIROSpectrum.h"
#include "KneiskeIROSpectrum.h"
#include "Stecker98IROSpectrum.h"
#include "Stecker2005IROSpectrum.h"
#include "Kneiske1001IROSpectrum.h"
#include "Kneiske0309IROSpectrum.h"
#include "Franceschini08EBL.h"
#include "Inoue12IROSpectrum.h"
#include "Medium.h"
#include "Stecker16Background.h"
#include "Franceschini17EBL.h"
#include "CustomBackground.h"

#define MAX_1D_COEF_NO 64
#define MAX_2D_COEF_NO 64

Medium::Medium():
neutrinoConc(0),
protonsConc(0),
m_pBackground(new CBackgroundTable(GetEBL())),
iBackgroundIntegral(0),
fB0(0),
fCurrentB(0),
fBincr(false),
neutrinoConc0(0),
fCoefTestEnabled(false)
{
	CouplingParameters params;
	fCoefTestEnabled = params.CoeffTestEnabled();
	double nuT=1.94858*units.K;// neutrino temperature at z=0 in internal units

	//neutrinoConc0=3/4/Pi/Pi*Zeta(3)*T^3
	neutrinoConc0=0.091345*nuT*nuT*nuT; //neutrino concentration at z=0 in internal units

	Reader()->readDoublePar("B0",fB0);
	fB0/=units.Bunit;
	Reader()->readBoolPar("B_INCREASING", fBincr);
	Reader()->readDoublePar("ProtonsConc",protonsConc);//read in cm^-3
	protonsConc *= units.Vunit;

	iBackgroundIntegral = new CBackgroundIntegral(*m_pBackground);
}

Medium::~Medium()
{
	delete iBackgroundIntegral;
}

int Medium::Update()
{
	m_pBackground->update(redshift.z());//this will also call iBackgroundIntegral->update(m_background);

	//background neutrino concentration
	neutrinoConc = neutrinoConc0 * redshift.dv();

	if(fCoefTestEnabled)
		m_pBackground->saveBackground("background");

	fCurrentB = fBincr?fB0*sqrt(redshift.dv()):fB0;

	return(0);
}

// IR/O background model to use
enum TIROBackgroundType
{
	IRO_None = 0,
	IRO_RESERVED1,
	IRO_RESERVED2,
	IRO_PRIMACK,
	IRO_KNEISKE,
	IRO_STECKER98,
	IRO_STECKER2005,
	IRO_BLACKBODY,
	IRO_KNEISKE_1001,//left for backward compatibility, use IRO_KNEISKE_MINIMAL (more precise definition)
	IRO_KNEISKE_0309,//left for backward compatibility, use IRO_KNEISKE_BEST_FIT (more precise definition)
	IRO_SARKAR_1005,
	IRO_KNEISKE_BEST_FIT,
	IRO_KNEISKE_MINIMAL,
	IRO_STECKER2012U,
	IRO_STECKER2012L,
	IRO_Inoue12_Baseline,
	IRO_Inoue12_LowerPop3,
	IRO_Inoue12_UpperPop3,
	IRO_Franceschini08,
	IRO_STECKER16Lower,
	IRO_STECKER16Upper,
    IRO_Franceschini17,

	IRO_END //must be the last
};

IBackgroundSpectrum* Medium::AddIRBebl()
{
	TIROBackgroundType IRO_model = IRO_STECKER2005;
	READ_SWITCH_SETTING(IRO_model, IRO_END);

	SmartPtr<IBackgroundSpectrum> iro;
	switch(IRO_model) {
	case IRO_PRIMACK:
		iro = new CPrimackIROSpectrum();
		break;
	case IRO_KNEISKE:
		iro = new CKneiskeIROSpectrum();
		break;
	case IRO_STECKER98:
		iro = new CStecker98IROSpectrum();
		break;
	case IRO_STECKER2005:
		iro = new CStecker2005IROSpectrum();
		break;
	case IRO_BLACKBODY:
	{
		double IRO_BB_T=11600;//Modified black body spectrum temperature (has effect if IRO_model==IRO_BLACKBODY)
		double IRO_BB_Power=0.;//Modified black body spectrum power (has effect if IRO_model==IRO_BLACKBODY)
		double IRO_BB_Density=-1.;//Modified black body background density in cm^-3 at z=0(has effect if IRO_model==IRO_BLACKBODY, if <=0 no normalizing is performed)
		bool   IRO_BB_Zdependence = false;//Modified black body background temperature Z-dependence if true T(Z) = (1+z)T(0)
		READ_DOUBLE_SETTING(IRO_BB_T);
		READ_DOUBLE_SETTING(IRO_BB_Power);
		READ_DOUBLE_SETTING(IRO_BB_Density);
		READ_BOOL_SETTING(IRO_BB_Zdependence);
		iro = new CModifiedBlackbodySpectrum(IRO_BB_T, IRO_BB_Power, IRO_BB_Density, IRO_BB_Zdependence, 1e100);
		break;
	}
	case IRO_KNEISKE_1001:
		iro = new Kneiske1001IROSpectrum();
		break;
	case IRO_KNEISKE_0309:
		iro = new Kneiske0309IROSpectrum();
		break;
	case IRO_SARKAR_1005:
		iro = new Sarkar1005Background(SarkarEBL);
		break;
	case IRO_KNEISKE_BEST_FIT:
		iro = new ElmagKneiskeBestFit();
		break;
	case IRO_KNEISKE_MINIMAL:
		iro = new ElmagKneiskeMinimal();
		break;
	case IRO_STECKER2012U:
		iro = new Stecker2012IROSpectrum(true);
		break;
	case IRO_STECKER2012L:
		iro = new Stecker2012IROSpectrum(false);
		break;
	case IRO_Inoue12_Baseline:
		iro = new Inoue12BaselineIROSpectrum();
		break;
	case IRO_Inoue12_LowerPop3:
		iro = new Inoue12LowPop3IROSpectrum();
		break;
	case IRO_Inoue12_UpperPop3:
		iro = new Inoue12UpperPop3IROSpectrum();
		break;
	case IRO_Franceschini08:
		iro = new Franceschini08EBL();
		break;
	case IRO_STECKER16Lower:
		iro = new Stecker16LowerBackground();
		break;
	case IRO_STECKER16Upper:
		iro = new Stecker16UpperBackground();
		break;
    case IRO_Franceschini17:
        iro = new Franceschini17EBL();
        break;
	case IRO_None:
		break;
	default:
		ThrowError("IR/O model chosen is not supported anymore");
	}
	if(iro!=0)
	{
		double iroMult = 1.;
		READ_DOUBLE_SETTING(iroMult);
		if(iroMult>0)
		{
			double IRO_extensionDeltaZconst = -1;//see HighRedshiftBackgrExtension
			READ_DOUBLE_SETTING(IRO_extensionDeltaZconst);
			if(IRO_extensionDeltaZconst>=0)
			{
				double IRO_extensionPowerLow = -3;//see HighRedshiftBackgrExtension
				double IRO_extensionDeltaZexp = 0;//see HighRedshiftBackgrExtension
				READ_DOUBLE_SETTING(IRO_extensionPowerLow);
				READ_DOUBLE_SETTING(IRO_extensionDeltaZexp);
				iro = new HighRedshiftBackgrExtension(iro, IRO_extensionDeltaZconst, IRO_extensionPowerLow, IRO_extensionDeltaZexp);
			}
			fEBL->addComponent(iro, iroMult);
		}
	}
    return iro;
}

// radio background model to use
enum TRadioBackgroundType
{
	Radio_None=0,
	Radio_Minimal,
	Radio_Minimal_Lee,//not supported anymore
	Radio_Middle,
	Radio_Maximal,
	Radio_Galactic,
	Radio_Sarkar1005,

	Radio_End  //must be the last
};

enum TCustomBackgroundUsage
{
    CustomBackgroundOff=0,
    CustomBackgroundAdd,
    CustomBackgroundReplace,

    CustomBackground_End
};

enum TCustomBackgroundEvolution
{
    CustomBackgroundEvolutionConstComoving=0,
    CustomBackgroundEvolutionProportional,
    CustomBackgroundEvolutionConstPhysical,

    CustomBackgroundEvolution_End
};

enum TCustomBackgroundFormat
{
    CustomBackgroundFormatPlain = 0,
    CustomBackgroundFormatMatrix,
    CustomBackgroundFormatAnalytic,

    CustomBackgroundFormatEnd
};

void Medium::AddRadioEbl()
{
	TRadioBackgroundType Radio_model = Radio_Middle;
	READ_SWITCH_SETTING(Radio_model, Radio_End);
	SmartPtr<IBackgroundSpectrum> radio;
	switch(Radio_model)
	{
	case Radio_None:
		break;
	case Radio_Minimal:
	case Radio_Galactic:
		radio = new DulkRadioBackground(Radio_model == Radio_Galactic);
		break;
	case Radio_Middle:
	case Radio_Maximal:
		radio = new ProtheroeRadioBackground(Radio_model == Radio_Maximal);
		break;
	case Radio_Sarkar1005:
		radio = new Sarkar1005Background(SarkarRadio);
		break;
	case Radio_Minimal_Lee:
		ThrowError("Minimal_Lee background support was dropped");
		break;
	default:
		ThrowError("Unexpected Radio_model parameter");
	}
	if(radio!=0)
	{
		double radioMult = 1.;
		READ_DOUBLE_SETTING(radioMult);
		if(radioMult>0)
			fEBL->addComponent(radio, radioMult);
	}
}

IBackgroundSpectrum* Medium::GetEBL()
{
	if(fEBL==0)
	{
        double BackgroundLowerCut_eV = 0.;
        double BackgroundUpperCut_eV = 1e100;
        READ_DOUBLE_SETTING(BackgroundLowerCut_eV);
        READ_DOUBLE_SETTING(BackgroundUpperCut_eV);
		fEBL = new CompoundBackground(BackgroundLowerCut_eV, BackgroundUpperCut_eV);
		bool NO_MWB = false;
		READ_BOOL_SETTING(NO_MWB);
		if(!NO_MWB)
		{
			double T0_CMB = 2.73;
			READ_DOUBLE_SETTING(T0_CMB);
			fEBL->addComponent(new CModifiedBlackbodySpectrum(T0_CMB, 0., -1., true, 1e100));
		}
		AddRadioEbl();
        AddIRBebl();
        TCustomBackgroundUsage customBackgroundUse = CustomBackgroundOff;

        if(Reader()->parExists("AddCustomBackgr", boolE)) //backward compatibility for settings file
        {
            bool AddCustomBackgr = false;
            READ_BOOL_SETTING(AddCustomBackgr);
            if(AddCustomBackgr)
                customBackgroundUse = CustomBackgroundAdd;
        }
        else {
            Reader()->readSwitchPar("CustomBackground", &customBackgroundUse, CustomBackground_End);
        }

		if(customBackgroundUse)
		{
            TCustomBackgroundEvolution CustomBackgroundEvolution = CustomBackgroundEvolutionConstComoving;
            READ_SWITCH_SETTING(CustomBackgroundEvolution, CustomBackgroundEvolution_End);
            TCustomBackgroundFormat CustomBackgroundFormat = CustomBackgroundFormatPlain;
            READ_SWITCH_SETTING(CustomBackgroundFormat, CustomBackgroundFormatEnd);

            CompoundBackground* evolution = 0;
            bool comoving = true;

            if(CustomBackgroundEvolution == CustomBackgroundEvolutionProportional)
                evolution = new CompoundBackground(fEBL);
            else if(CustomBackgroundEvolution == CustomBackgroundEvolutionConstPhysical)
                comoving = false;

//char custom_backgr_file[DEFAULT_STR_PAR_MAX_LENGTH] = "";
            std::string custom_backgr_file = "tables/custom_ebl";
			READ_STRING_SETTING(custom_backgr_file);
			SmartPtr<IBackgroundSpectrum> b = 0;
			if(CustomBackgroundFormat == CustomBackgroundFormatMatrix) {
				b = new MatrixBackground(custom_backgr_file, comoving, false);
			}
			else if(CustomBackgroundFormat == CustomBackgroundFormatPlain) {
                /// if predefined evolution is set always use comoving mode for selfconsistency
				b = new PlainTableBackground(custom_backgr_file, comoving || evolution, evolution);
			}
            else if (CustomBackgroundFormat == CustomBackgroundFormatAnalytic){
                b = new CustomBackground();
            }
            else
                ThrowError("unsupported CustomBackgroundFormat");

            if(customBackgroundUse == CustomBackgroundAdd)
			    fEBL->addComponent(b);
            else if(customBackgroundUse == CustomBackgroundReplace)
                fEBL->replacePart(b);
            else
                NOT_IMPLEMENTED;

			CopyFileOrDir(custom_backgr_file, plt_local_dir + DIR_DELIMITER_STR);
		}
		fEBL->init();
        double density = BackgroundUtils::CalcIntegralDensity(*fEBL, 0);
        cerr << "background density at z=0 [cm^-3]: " << density << endl;
	}
	return fEBL;
}

SmartPtr<CompoundBackground> Medium::fEBL;

double Medium::neutrinoClusteringModifier=1.0;//current clustering modifier (=1. no clustering)

