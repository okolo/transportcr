#ifndef StringResource_H_INCLUDED
#define StringResource_H_INCLUDED

#include <string>

class StringResource
{
public:
	//construct empty resource
	StringResource();
	//constructur with base64 encoded string array
	StringResource(const std::string* aStringArray, int aDataLength);
	std::string GetString();
	inline int IsEmpty() {return fStringArray == 0;};
private:
	const std::string* fStringArray;
	int fDataLength;
};

#endif //#ifndef StringResource_H_INCLUDED
