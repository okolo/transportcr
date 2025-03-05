#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "Addfunc.h"
#include <stdio.h>
#include <signal.h>

//#include "Prodspec.h"
#include <string.h>
#include "time.h"

#include "main.h"
#include "VersionInfo.h"
#include "const.h"
#include "Concentrations.h"
#include "TimeZ.h"
#include "LEPData.h"
#include "KneiskeIROSpectrum.h"
#include "XMLSwitchFileReader.h"
#include "Units.h"
#include "Ranges.h"
#include "PropagEngine.h"
#include "ParticleList.h"
#include "PowerLowInjectionSpectra.h"
#include "TableInjSpectra.h"
#include "CustomInjSpectrum.h"
#include "FragmentationBasedInjSpectra.h"
#include "BlackHoleInjSpectra.h"
#include "Background.h"
#include "ScanInfoWritter.h"
#include "Jacobian.h"
#include "Coupling.h"
#include "MassiveNeutrino.h"
#include "Galaxy.h"
#include "BLLac.h"

void SIGSEGV_handler(int sig) {
    print_stack_trace();
    exit(1);
}

void Application::GslFrrorHandler (const char * reason,
        const char * file,
        int line,
        int gsl_errno)
{
	std::ostringstream message;
	message << "GSL error " << gsl_errno << " occurred in " << file << "("<< line << ")\nreason:" << reason << "\n";
	std::string msg = message.str();
	ASSERT(false);
	ThrowError(msg);
}

Application::Application()
{
	gsl_set_error_handler (GslFrrorHandler);
}

const char* Application::progName=0;

void Application::usageExit()
{
	cerr << "Usage:"<< endl;

	cerr << progName << " file.xsw" << endl
		<< "-run using parameter file in XML format file.xsw" << endl <<endl
		<< progName << " --printPatch file.diff" << endl
		<< "print git patch file (useful for dirty version builds)" << endl;
	Exit(retUsageError);
}


/*88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
88888888888888888888888888888888888888     M A I N    88888888888888888888888888888888888888*/

string plt_local_dir(PLT_DIR_LINK DIR_DELIMITER_STR "tmp_plt");//local directory to store temporary results
const char* plt_local_c;

int main(int argc,char *argv[])
{
    signal(SIGSEGV, SIGSEGV_handler);
	return Application::main(argc, argv);
}

int Application::main(int argc,char *argv[])
{
	VersionInfo::PrintVersionInfo(cerr, "PROPAGATION");

	progName = GetFileName(argv[0]);
	const char* xmlExt = ".xsw";
	
	if(argc==1)
		usageExit();


	if(strcmp(argv[1],"--printPatch")==0)
	{
		if(argc!=3)
			usageExit();
		VersionInfo::PrintDiffFile(argv[2]);
		return 0;
	}
	if(strcmp(argv[1],"--test")==0)
	{
		return UnitTest(argc-2, argv + 2);
	}

	string parFile = GetFileName(argv[1]);
	int pointPos = parFile.find_last_of('.');
	if (pointPos<1) {
		usageExit();
	}
	string parFileName(parFile.c_str(),pointPos);
	string parFileExt(parFile.c_str()+pointPos);
	string pltDir(ALL_PLT_DIR);
	pltDir += parFileName;
	plt_local_dir += parFileName;
	plt_local_dir += DIR_DELIMITER_STR;
	plt_local_c = plt_local_dir.c_str();
	string switchFileCopy = plt_local_dir + parFile;
	RemoveFileOrDir(plt_local_c);
	Mkdir(plt_local_c);
	VersionInfo::PrintDiffFile((plt_local_dir + "propag.diff").c_str());
	try{
		if(parFileExt == xmlExt || parFileExt == ".xml")
		{
			LOG_MESSAGE2("Processing ", argv[1]);
			const char* dtdFile = "switches.dtd";
			CopyFileOrDir(dtdFile, plt_local_dir + dtdFile);
			CXMLSwitchFileReader xmlReader(argv[1]);
			if(!xmlReader.parse())
				Exit(retInvalidSettingsFile);
			time_t t = time(0);
			Parameters::SetReader(&xmlReader);
			IParWriter::Set(&xmlReader);//may be used to output adjusted param values
			VersionInfo::SaveVersionInfo(&xmlReader);
			Application app;
			app.run();
			xmlReader.print();//saving version info and adding omitted parameters
			t = (time(0)-t);
			LOG_MESSAGE3("Calculation took ", t, " sec");
		}
		else
			usageExit();
	}catch (const char* errMsg) {
		cerr << errMsg << endl;
        print_saved_stack_trace();
	}

	CopyFileOrDir(argv[1],switchFileCopy);
	RemoveFileOrDir(pltDir);
	if(!MoveFileOrDir(plt_local_c,pltDir))
	{
		string erMsg = "failed to move " + plt_local_dir + " to " + pltDir;
		ReportError(erMsg.c_str());
	}

	return 0;
}

int Application::UnitTest(int argc,char *argv[])
{
	try{
	if(!argc)
	{
		//TODO: default test
		return 0;
	}
	std::string testName = argv[0];
	if(testName=="BLLac")
	{
		BLLacLDDE::UnitTest();
	}
	}catch(const char* errMsg)
	{
		std::cerr << errMsg << std::endl;
		return 1;
	}
	return 0;
}

CInjectionSpectra* Application::CreateSource()
{
	enum TInjectionSpectra{
		EPowerLaw=0,
		EFixedEnergy,
		ECustomTable,
		EFragmentationBased,
		ENeronovBH,
        ECustomCode,

		ENumberOfInjSpecImpl
	};
	TInjectionSpectra specType = EPowerLaw;
	Reader()->readSwitchPar("InjectionSpectraType",&specType,ENumberOfInjSpecImpl);
	switch(specType)
	{
	case EPowerLaw:
		return new CPowerLowInjectionSpectra();
	case EFixedEnergy:
		return new CFixedEnergyInjection();
	case ECustomTable:
		return new CTableInjSpectra();
	case EFragmentationBased:
		return new CFragmentationBasedInjSpectra();
	case ENeronovBH:
		return new CBlackHoleInjSpectra();
	case ECustomCode:
	    return new CCustomInjSpectrum();
	default:
		ThrowError("Unexpected InjectionSpectraType");
	}
	return 0;
}

enum EngineMode
{
	LongDistance = 0,
	ShortDistance,
	GalaxyInnerSpec,
	GalaxyOuterSpec,
	End_EngineMode
};

void Application::run()
{
	EngineMode mode = LongDistance;
	double MinInteractionEnergyMeV = 1.0;//MeV  the energy below which all particles do not interact
    double MinInteractionEnergyNucleiMeV = 0; //MeV the energy per nucleus below which nuclei do not interact
	int energyResolutionFactor = 2;
	double EminEv = 1e6;
	Reader()->readDoublePar("Emin",EminEv);
	READ_INT_SETTING(energyResolutionFactor);

	int zResolutionFactor = 1;
	READ_INT_SETTING(zResolutionFactor);

	Reader()->readSwitchPar("CalcMode",&mode,End_EngineMode);

	READ_DOUBLE_SETTING(MinInteractionEnergyMeV);
    READ_DOUBLE_SETTING(MinInteractionEnergyNucleiMeV);
	double H_in_km_s_Mpc = 71;
	READ_DOUBLE_SETTING(H_in_km_s_Mpc);
	double Lv = 0.73;//default value
	READ_DOUBLE_SETTING(Lv);

	SmartPtr<CInjectionSpectra> spec = CreateSource();
	double maxZ = spec->MaxZ();
	double Emax = spec->BaseEmax()/units.Eunit;

	CRanges::SetDefault(
		new CRanges(
			energyResolutionFactor,
			EminEv*1e-6/units.Eunit,
			Emax,
			maxZ,
			zResolutionFactor,
			CPropagEngine::useCELredshift(),
			MinInteractionEnergyMeV,
            MinInteractionEnergyNucleiMeV
			)
		);

	CTimeZ::init(Lv,H_in_km_s_Mpc);
	redshift.setZ(0);

	spec->init();
	spec->NormalizeIfNeeded();
	CouplingList::Instance()->init();

	redshift.setZ(Ranges().Zmax());
	std::cout << "Zmax=" << Ranges().Zmax() << " (L=" << CTimeZ::z2d(Ranges().Zmax())*units.Lunit/units.Mpc_cm << " Mpc)\n";

	double microStep=0.5;//here in Mpc, will be converted to internal units
	READ_DOUBLE_SETTING(microStep);
	microStep*=units.Mpc_cm /units.Lunit;
	switch(mode)
	{
	case LongDistance:
		{
			Concentrations N;
			CLongDistanceEngine	engine(Ranges().Z(), spec);
			TPropagCoef c = engine.runZ(N, microStep);
			c.free();
			N.PrintSpectrum("uniform");
		}
		break;
	case ShortDistance:
		{
			redshift.setZ(0.);
			Concentrations N1;
			CShortDistanceTestEngine shortDistanceTestEngine(spec);
			shortDistanceTestEngine.run(N1, microStep);
		}
		break;
	case GalaxyInnerSpec:
		{
			double GalL = 0.;
			double GalB = -90;
			double GalDeltaL = 360.;
			double GalDeltaB = 180.;
			READ_DOUBLE_SETTING(GalL);
			READ_DOUBLE_SETTING(GalB);
			READ_DOUBLE_SETTING(GalDeltaL);
			READ_DOUBLE_SETTING(GalDeltaB);
			Concentrations N;
			GalaxyPropagEngine	engine(spec, microStep);
			engine.run(N, GalL, GalL + GalDeltaL, GalB, GalB + GalDeltaB);
			N.PrintSpectrum("galaxy");
		}
		break;
	case GalaxyOuterSpec:
		{
			Concentrations N;
			GalaxyEffectiveSourceEngine	engine(spec, microStep);
			engine.run(N);
			N.PrintSpectrum("galaxy_spec");
		}
		break;
	default:
		ThrowError("Invalid CalcMode parameter");
	}
};
