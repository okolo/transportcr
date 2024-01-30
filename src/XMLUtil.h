// XMLUtil.h: interface for the CXMLUtil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_XMLUTIL_H__F638E4A8_0348_4885_AFE1_F5DC7EFEF6A1__INCLUDED_)
#define AFX_XMLUTIL_H__F638E4A8_0348_4885_AFE1_F5DC7EFEF6A1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//using namespace std;

#include <iostream>
#include <iosfwd>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/util/TranscodingException.hpp>


//#include <xercesc/dom/deprecated/DOM_DOMException.hpp>
//#include <xercesc/dom/deprecated/DOMStringWrapper.hpp>

//#include <parsers/DOMParser.hpp>
//#include <xercesc/dom/deprecated/DOMParser.hpp>
//#include <xercesc/dom/deprecated/DOM.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMNode.hpp>

XERCES_CPP_NAMESPACE_USE

class DOMStringWrapper
{
public:
	DOMStringWrapper(const XMLCh* aString);
	DOMStringWrapper(const char* aString);
	DOMStringWrapper(XMLCh* aString);
	DOMStringWrapper(char* aString);
	~DOMStringWrapper();
	operator const char* () const;
	operator const XMLCh*() const;
	bool operator==(const XMLCh*) const;
	bool operator==(const char*) const;
	bool operator==(const DOMStringWrapper&) const;
	int length() const;
private:
	XMLCh* lString;
	char* sString;
	bool fOwnsL;
	bool fOwnsS;
};

///The class is used for automatic xerces lib (un)initialization
class CXMLUtil  
{
public:
	CXMLUtil();
	virtual ~CXMLUtil();
private:
	static int s_refNo;
};

//ostream& operator<< (ostream& target, const DOMStringWrapper& s);
//XMLFormatter& operator<< (XMLFormatter& strm, const DOMStringWrapper& s);
//ostream& operator<<(ostream& target, DOMNode& toWrite);


#endif // !defined(AFX_XMLUTIL_H__F638E4A8_0348_4885_AFE1_F5DC7EFEF6A1__INCLUDED_)
