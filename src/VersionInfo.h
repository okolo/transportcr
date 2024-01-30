#ifndef VersionInfo_INCLUDED
#define VersionInfo_INCLUDED

#include <iostream>
#include <string>

class IParWriter;

class VersionInfo
{
public:
#ifdef USE_SubWCRev
    static const char* WCREV;
    static const char* WCDATE; 
    static const char* WCNOW;
    static const char* WCRANGE;
    static const char* WCMIXED;
    static const char* WCMODS;
    static const char* WCURL;
#endif
	static void PrintDiffFile(const char* FilePath);
	static void PrintDiff(std::ostream& aOut); 
	static void PrintVersionInfo(std::ostream& aOut, std::string aProgramName);
	static void SaveVersionInfo(IParWriter* aWriter);
};

#endif //#ifndef VersionInfo_INCLUDED
