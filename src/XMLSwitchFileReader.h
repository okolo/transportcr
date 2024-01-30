// XMLSwitchFileReader.h: interface for the CXMLSwitchFileReader class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_XMLSWITCHFILEREADER_H__BF210CD1_EE84_11D5_B7D3_93CED5480B9A__INCLUDED_)
#define AFX_XMLSWITCHFILEREADER_H__BF210CD1_EE84_11D5_B7D3_93CED5480B9A__INCLUDED_

#include "XMLUtil.h"

#include <string.h>
#include <stdlib.h>
#include "Parameters.h"
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

class CXMLSwitchFileReader  :  public CXMLUtil , public IParReader, public IParWriter
{
public:
	CXMLSwitchFileReader(const char* theXmlFile);
	virtual ~CXMLSwitchFileReader();

//	implemented from IParReader

	bool parse();

    bool parExists(const char* theParName, EParameterType theParType);

	bool readDoublePar( const char* theParName, double& thePar);

	bool readIntPar( const char* theParName, int& thePar);

	bool readBoolPar( const char* theParName, bool& thePar);

	bool readSwitchPar( const char* theParName, int& thePar);
				   
	bool readStringPar( const char* theParName, char* thePar, int maxLength = DEFAULT_STR_PAR_MAX_LENGTH);

	inline void redirectToStdOut(bool aRedirect = true) {m_redirectedOutput = aRedirect;}

	bool save(const char* theXmlFile);

	// implemented from IParWriter
	bool print();

	bool writePar( const char* theParName, double thePar, const char *theComment=0);

	bool writePar( const char* theParName, int thePar, const char *theComment=0);

	bool writeBoolPar( const char* theParName, bool thePar, const char *theComment=0);

	bool writeSwitchPar( const char* theParName, int thePar, const char *theComment=0);

	bool writePar( const char* theParName, const char* thePar, const char *theComment=0);

	void WriteLine(const char* theComment = 0){};

	bool SupportsPositioning() {return true;};

protected:

	bool createPar(const char* aParName,
		const char* aParValue,
		EParameterType aParType, const char* aComment);


	bool writePar(const char* aParName,
		const char* aParValue,
		EParameterType aParType, const char* aComment);

	//Throws theParName if can't find the parameter
	DOMNode* findPar(const char* aParName, EParameterType aParType, const char* aDefaultValue);

	//Throws theParName if can't find the parameter
	DOMElement* findParTag(const char* aParName, EParameterType aParType);
/**
Throws theParName if can't find the parameter
the returned string must be deleted
*/
	const XMLCh* readPar(const char* aParName, EParameterType aParType, std::string aDefaultValue);
	void addMissingPar(const char* aParName, EParameterType aParType, std::string aValue, const char* aComment=0);
	
	bool readIntegerPar( const char* aParName, int& aPar, EParameterType aParType);

	bool serializeDOM(DOMNode* node, const char* aFileName);

	const char* convert(int aVal);
	const char* convert(double aVal);
	const char* convertbool(bool aVal);
	
	ErrorHandler *m_erHandler;
	const char* m_XmlFile;
	XercesDOMParser *m_parser;
	XERCES_CPP_NAMESPACE::DOMDocument* m_doc;
	DOMNodeList* m_parArrays[endParameterTypeE];
	static const char* sTypeNames[endParameterTypeE];//the array of parameter names
	DOMNode*  m_newParDir;//used for writing previously undefined settings
	static const char* sParNameAttr;
	static const char* sParValueAttr;
	static const char* sParCommentAttr;
	static const char* sParDefaultValueAttr;
	char				m_buffer[64];
	bool				m_redirectedOutput;
};

#endif // !defined(AFX_XMLSWITCHFILEREADER_H__BF210CD1_EE84_11D5_B7D3_93CED5480B9A__INCLUDED_)
