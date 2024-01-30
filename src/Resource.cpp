#include "Resource.h"
#include "base64.h"

StringResource::StringResource(const std::string* aStringArray, int aDataLength):
fStringArray(aStringArray),
fDataLength(aDataLength)
{
}

StringResource::StringResource():
fStringArray(0),
fDataLength(0)
{
}

std::string StringResource::GetString()
{
	std::string tot = "";
	for(int i=0; i<fDataLength; i++)
		tot = tot + base64_decode(fStringArray[i]);
	return tot;
}
