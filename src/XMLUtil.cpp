// XMLUtil.cpp: implementation of the CXMLUtil class.
//
//////////////////////////////////////////////////////////////////////

#include "XMLUtil.h"
#include <iostream>
using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int CXMLUtil::s_refNo = 0;

CXMLUtil::CXMLUtil()
{
	if(s_refNo==0)
	{
		try
	    {
	        XMLPlatformUtils::Initialize();
	    }
	    catch(const XMLException& toCatch)
	    {
	        cerr << "Error during Xerces-c Initialization.\n"
	             << "  Exception message:"
	             << (const char*)DOMStringWrapper(toCatch.getMessage()) << endl;
	        throw &toCatch;
	    }
	}
	s_refNo++;
}

CXMLUtil::~CXMLUtil()
{
	s_refNo--;
	if(s_refNo==0)
	{
		XMLPlatformUtils::Terminate();
	}
}

DOMStringWrapper::DOMStringWrapper(const XMLCh* aString):
	lString((XMLCh*)aString),
	sString(0),
	fOwnsL(0),
	fOwnsS(1)
{
	sString = XMLString::transcode(aString);
}

DOMStringWrapper::DOMStringWrapper(const char* aString):
	lString(0),
	sString((char*)aString),
	fOwnsL(1),
	fOwnsS(0)
{
	lString = XMLString::transcode(aString);
}

DOMStringWrapper::DOMStringWrapper(XMLCh* aString):
	lString(aString),
	sString(0),
	fOwnsL(1),
	fOwnsS(1)
{
	sString = XMLString::transcode(aString);
}

DOMStringWrapper::DOMStringWrapper(char* aString):
	lString(0),
	sString(aString),
	fOwnsL(1),
	fOwnsS(1)
{
	lString = XMLString::transcode(aString);
}

DOMStringWrapper::operator const char* () const
{
	return sString;
}

DOMStringWrapper::operator const XMLCh*() const
{
	return lString;
}

int DOMStringWrapper::length() const
{
	return strlen(sString);
}

DOMStringWrapper::~DOMStringWrapper()
{
	if(fOwnsL && lString)
		XMLString::release(&lString);
	if(fOwnsS && sString)
		XMLString::release(&sString);
}

bool DOMStringWrapper::operator==(const XMLCh* aStr) const
{
	return XMLString::equals(lString, aStr);
}

bool DOMStringWrapper::operator==(const char* aStr) const
{
	return XMLString::equals(sString, aStr);
}

bool DOMStringWrapper::operator==(const DOMStringWrapper& aStr) const
{
	return XMLString::equals(lString, aStr.lString);
}
