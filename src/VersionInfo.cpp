#include "Addfunc.h"
#include "VersionInfo.h"
#include "Resource.h"
#include <fstream>
#include "Parameters.h"

#include "_autogeneratedVersionInfo.h"

void VersionInfo::PrintVersionInfo(std::ostream& aOut, std::string aProgramName)
{
	aOut << "\n" << aProgramName << " " << RELEASE_INFO << "\n" <<
#ifdef _DEBUG
	"debug"
#else
	"release"
#endif
	<< " version" << std::endl;
}

#include "_autogeneratedGitInfo.c"
#include "_autogeneratedDif.c"

void VersionInfo::SaveVersionInfo(IParWriter* aWriter)
{
	extern StringResource EncodedVersion;
	extern StringResource EncodedDif;
	std::string ver = "0";
	int difHash = -1;
	if(!EncodedVersion.IsEmpty())
	{
		difHash = 0;
		ver = EncodedVersion.GetString();
		if(!EncodedDif.IsEmpty())
		{
			std::string decoded = EncodedDif.GetString();
			difHash = (int)hash_string(decoded.c_str());
		}
	}
	aWriter->writePar("Revision", ver.c_str(), "version info obtained with svnversion command or 0 if it is not available");
	aWriter->writePar("DifHash", difHash, "diff file hash or 0 if it is clean version or -1 if version info is not available");
}

void VersionInfo::PrintDiff(std::ostream& aOut)
{
	extern StringResource EncodedDif;
	if(EncodedDif.IsEmpty())
		return;
	std::string decoded = EncodedDif.GetString();
	aOut << decoded;
}

void VersionInfo::PrintDiffFile(const char* FilePath)
{
	extern StringResource EncodedDif;
	if(EncodedDif.IsEmpty())
		return;//don't create empty file if no local changes found

	std::ofstream file(FilePath, std::ios_base::out | std::ios_base::binary);
	if(!file.is_open())
		throw "PrintDiffFile: failed to open file for writing";
	PrintDiff(file);
	file.close();
}

