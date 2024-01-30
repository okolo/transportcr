// XMLSwitchFileReader.cpp: implementation of the CXMLSwitchFileReader class.
//
//////////////////////////////////////////////////////////////////////

#include "XMLSwitchFileReader.h"
#include <fstream>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include "Addfunc.h"

const char* CXMLSwitchFileReader::sTypeNames[endParameterTypeE] = {
//order of elements the array should be synchronized with EParameterType values
		"double",//doubleE
		"int",//intE
		"bool",//boolE
		"switch",//switchE
		"string"//stringE
};

const char* CXMLSwitchFileReader::sParNameAttr = "title";
const char* CXMLSwitchFileReader::sParValueAttr = "value";
const char* CXMLSwitchFileReader::sParCommentAttr = "comment";
const char* CXMLSwitchFileReader::sParDefaultValueAttr = "default";

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
//  ostream << DOMStringWrapper
//
//  Stream out a DOM string. Doing this requires that we first transcode
//  to char * form in the default code page for the system
// ---------------------------------------------------------------------------

CXMLSwitchFileReader::CXMLSwitchFileReader(const char* theXmlFile):
m_XmlFile(theXmlFile),
m_redirectedOutput(false)
{
	/*
    m_parser = new DOMParser;
    m_parser->setValidationScheme(DOMParser::Val_Auto);
    m_parser->setDoNamespaces(false);
    m_parser->setDoSchema(false);
    m_parser->setValidationSchemaFullChecking(false);
    m_erHandler = new DOMTreeErrorReporter();
    m_parser->setErrorHandler(m_erHandler);
    m_parser->setCreateEntityReferenceNodes(false);
    m_parser->setToCreateXMLDeclTypeNode(true);
	*/
    m_parser = new XercesDOMParser();
    m_parser->setValidationScheme(XercesDOMParser::Val_Auto);
    m_parser->setDoNamespaces(false);    // optional

    m_erHandler = (ErrorHandler*) new HandlerBase();
    m_parser->setErrorHandler(m_erHandler);
}

CXMLSwitchFileReader::~CXMLSwitchFileReader()
{
	delete m_parser;
	delete m_erHandler;
}

bool CXMLSwitchFileReader::parse()
{
  //bool retval = true;
	bool errorsOccured = true;
    try {
        m_parser->parse(m_XmlFile);
		errorsOccured = false;
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cerr << "XML parse error: \n"
             << message << "\n";
        XMLString::release(&message);
    }
    catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cerr << "XML parse error: \n"
             << message << "\n";
        XMLString::release(&message);
    }
    catch (const SAXException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cerr << "XML parse error: \n"
             << message << "\n";
        XMLString::release(&message);
    }
    catch (...) {
        cerr << "Unknown XML parse error occurred \n" ;
    }

    // If the parse was successful, output the document data from the DOM tree
    if (errorsOccured)
		return false;

    m_doc = m_parser->getDocument();

	for(int i=0; i<endParameterTypeE; i++)
	  {
	    DOMStringWrapper typeName(sTypeNames[i]);
	    m_parArrays[i] = m_doc->getElementsByTagName(typeName);
	  }

	m_newParDir = m_doc->getFirstChild();
    return true;
}

/**
Throws theParName if can't find the parameter
*/
const XMLCh* CXMLSwitchFileReader::readPar(const char* aParName, 
		EParameterType aParType, std::string aDefaultValue)//using std::string to ensure that local copy of argument is used
{
	const XMLCh* value = NULL;
	try{
		DOMNode* node = findPar(aParName, aParType, aDefaultValue.c_str());
		value = node->getNodeValue();
	}catch(const char* str){
		addMissingPar(aParName, aParType, aDefaultValue);
	}
	return value;
}

void CXMLSwitchFileReader::addMissingPar(const char* aParName, EParameterType aParType, std::string aValue, const char* aComment)
{
	const char* MissingParamsCountTitle = "MissingParamsCount";
	if(strcmp(aParName,MissingParamsCountTitle)) {//avoid infinite loop
		int missingParCount = 0;
		readIntPar("MissingParamsCount", missingParCount);//this will position m_newParDir to the right place
		missingParCount++;
		writePar("MissingParamsCount", missingParCount, "number of missing parameters");
	}
	else {//if MissingParamsCountTitle tag doesn't exist
		try {
			findParTag("scanPar", doubleE);//try to position new tag next to "scanPar" in "advanced" section
		} catch (const char *) { }

	}
	std::cerr << "creating new entry for missing parameter " << aParName << "=" << aValue << std::endl;
	createPar(aParName, aValue.c_str(), aParType, aComment);
}

/* old variant
char* CXMLSwitchFileReader::readPar(const char* theParName, 
		EParameterType theParType)
{
	ASSERT(theParType!=endParameterTypeE);

	DOMStringWrapper parNameStr(theParName);
	DOMStringWrapper titleStr("title");
	DOMStringWrapper valStr("value");

	DOMNodeList list = m_parArrays[theParType];
	unsigned int iMax = list.getLength();
	for(unsigned int i = 0; i<iMax; i++)
	{
		DOMNode el = list.item(i);
		DOM_NamedNodeMap map = el.getAttributes();
		DOMStringWrapper parName = map.getNamedItem(titleStr).getNodeValue();
		if(parNameStr.equals(parName))
		{
			char* value = DOMString2str(map.getNamedItem(valStr).getNodeValue());//getAttribute(valStr));
			return value;
		}
	}
	throw theParName;
}
*/


bool CXMLSwitchFileReader::parExists(const char* theParName, EParameterType theParType){
	try{
		findParTag(theParName, theParType);
	}
	catch (const char*){
		return false;
	}
	return true;
}

/**
Find parameter and set it's default value (if aDefaultValue!=NULL)
Throws theParName if can't find the parameter
*/
DOMNode* CXMLSwitchFileReader::findPar(const char* theParName, EParameterType theParType, const char* aDefaultValue)
{
	ASSERT(theParType!=endParameterTypeE);
	DOMStringWrapper titleStr(sParNameAttr);
	DOMStringWrapper valStr(sParValueAttr);

	DOMElement* el = findParTag(theParName, theParType);
	if(aDefaultValue!=NULL)
	{
		DOMStringWrapper attr(sParDefaultValueAttr);
		DOMStringWrapper val(aDefaultValue);
		el->setAttribute(attr, val);
	}
	return el->getAttributes()->getNamedItem(valStr);
}

DOMElement* CXMLSwitchFileReader::findParTag(const char* theParName, EParameterType theParType)
{
	ASSERT(theParType!=endParameterTypeE);

	DOMStringWrapper parNameStr(theParName);
	DOMStringWrapper titleStr(sParNameAttr);
	DOMStringWrapper valStr(sParValueAttr);

	//	((DOMElement)node).setAttribute(defaultValAttr, defaultParValue);

	DOMNodeList* list = m_parArrays[theParType];
	unsigned int iMax = list->getLength();
	for(unsigned int i = 0; i<iMax; i++)
	{
		DOMElement* el = (DOMElement*)list->item(i);
		const XMLCh* parName = el->getAttribute (titleStr);

		if(parNameStr==parName)
		{
			m_newParDir = el->getParentNode();
			return el;
		}
	}
	throw theParName;
}

bool CXMLSwitchFileReader::readDoublePar( const char* theParName,
				   double& thePar)
{
	const XMLCh* val = readPar(theParName, doubleE, convert(thePar)); 
	if(val==NULL)
		return false;
	DOMStringWrapper value = val;
	bool result = true;
	char* endPtr;
	double par = strtod(value,&endPtr);
	if(endPtr[0]!='\0')
		result = false;
	else
		thePar = par;
	return result;
}

bool CXMLSwitchFileReader::readIntegerPar( const char* theParName,
				   int& thePar, EParameterType theParType)
{
	const XMLCh* val = readPar( theParName, theParType, convert(thePar)); 
	if(val==NULL)
		return false;
	DOMStringWrapper value = val;
	bool result = true;
	char* endPtr;
	unsigned long par = strtoul(value,&endPtr,10);
	if(endPtr[0]!='\0')
		result = false;
	else
		thePar = (int)par;
	return result;
}

bool CXMLSwitchFileReader::readIntPar( const char* theParName,
				   int& thePar)
{
	return readIntegerPar( theParName, thePar, intE);
}

bool CXMLSwitchFileReader::readBoolPar( const char* theParName, bool& thePar)
{
	const XMLCh* val = readPar( theParName, boolE, convertbool(thePar));

	if(val==NULL)
		return false;
	DOMStringWrapper value = val;
	bool result = true;

	if(value == "true")
		thePar = true;
	else if(value == "false")
		thePar = false;
	else
		result = false;

	return result;
}

bool CXMLSwitchFileReader::readSwitchPar( const char* theParName,
				   int& thePar)
{
	return readIntegerPar( theParName, thePar, switchE); 
}

bool CXMLSwitchFileReader::save(const char* theXmlFile)
{
	return serializeDOM(m_doc, theXmlFile);
}

bool CXMLSwitchFileReader::serializeDOM(DOMNode* node, const char* aFileName) 
{
    XMLCh tempStr[1000];
    XMLString::transcode("LS", tempStr, 999);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);
    DOMLSSerializer* theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();

    // optionally you can set some features on this serializer
    if (theSerializer->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
        theSerializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);

    if (theSerializer->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
         theSerializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);

    // optionally you can implement your DOMLSSerializerFilter (e.g. MyDOMLSSerializerFilter)
    // and set it to the serializer
    //DOMLSSerializer* myFilter = new myDOMLSSerializerFilter();
    //theSerializer->setFilter(myFilter);

    // optionally you can implement your DOMErrorHandler (e.g. MyDOMErrorHandler)
    // and set it to the serializer
    //DOMErrorHandler* errHandler = new myDOMErrorHandler();
    //theSerializer->getDomConfig()->setParameter(XMLUni::fgDOMErrorHandler, myErrorHandler);

    // StdOutFormatTarget prints the resultant XML stream
    // to stdout once it receives any thing from the serializer.
	XMLFormatTarget *myFormTarget = 0;
	DOMLSOutput* theOutput = 0;
	const char* errMsg = "Failed to save settings to xml: \n";
	try {
		myFormTarget = 
			(aFileName==0)?
				((XMLFormatTarget*)new StdOutFormatTarget()):
				((XMLFormatTarget*)new LocalFileFormatTarget(aFileName));
		theOutput = ((DOMImplementationLS*)impl)->createLSOutput();
		theOutput->setByteStream(myFormTarget);

        // do the serialization through DOMLSSerializer::write();
        theSerializer->write(node, theOutput);
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cerr << errMsg
             << message << "\n";
        XMLString::release(&message);
        return false;
    }
    catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cerr << errMsg
             << message << "\n";
        XMLString::release(&message);
        return false;
    }
    catch (...) {
        cerr << errMsg << "Unexpected Exception \n" ;
        return false;
    }

    theOutput->release();
    theSerializer->release();
    delete myFormTarget;
    return true;
}


bool CXMLSwitchFileReader::print()
{
	if(m_redirectedOutput)
	{
		return serializeDOM(m_doc, 0);
	}
	else
		return save(m_XmlFile);
}

bool CXMLSwitchFileReader::writePar(const char* theParName,
		const char* theParValue,
		EParameterType theParType, const char* aComment)
{
	try{
		DOMNode* node = findPar(theParName, theParType, NULL);
		DOMStringWrapper val(theParValue);
		node->setNodeValue(val);
	}
	catch(const char* name)
	{
		addMissingPar(theParName, theParType, theParValue, aComment);
		return false;
	}
	return true;
}

bool CXMLSwitchFileReader::createPar(const char* aParName,
		const char* aParValue,
		EParameterType aParType, const char* aComment)
{
	if (aComment==NULL) {
		aComment = aParName;
	};

	DOMStringWrapper titleAttr(sParNameAttr);
	DOMStringWrapper valueAttr(sParValueAttr);
	DOMStringWrapper commentAttr(sParCommentAttr);
	DOMStringWrapper defaultValAttr(sParDefaultValueAttr);

	DOMStringWrapper theParName(aParName);
	DOMStringWrapper theParValue(aParValue);
	DOMStringWrapper theParComment(aComment);
	DOMStringWrapper tagName(sTypeNames[aParType]);
	try{
		DOMElement* newPar = m_doc->createElement(tagName);
		newPar->setAttribute(titleAttr, theParName);
		newPar->setAttribute(valueAttr, theParValue);
		newPar->setAttribute(defaultValAttr, theParValue);
		newPar->setAttribute(commentAttr, theParComment);

		m_newParDir->appendChild(newPar);
		m_parArrays[aParType] = m_doc->getElementsByTagName(tagName);//make sure future searches ot the parameter will succeed
	}catch (const char* name) {
		char erMsg[256];
		sprintf(erMsg,"failed to create tag for \"%s\" setting",name);
			ReportError(erMsg);
		return false;
	}

	return true;
}

const char* CXMLSwitchFileReader::convert(int aVal)
{
	sprintf(m_buffer,"%d",aVal);
	return m_buffer;
}

const char* CXMLSwitchFileReader::convert(double aVal)
{
	sprintf(m_buffer,"%lg",aVal);
	return m_buffer;
}

const char* CXMLSwitchFileReader::convertbool(bool aVal)
{
	return (aVal)?"true":"false";
}

bool CXMLSwitchFileReader::writePar( const char* theParName,
	   double thePar, const char *theComment)
{
	return writePar(theParName, convert(thePar), doubleE, theComment);
}

bool CXMLSwitchFileReader::writePar( const char* theParName,
	   int thePar, const char *theComment)
{
	return writePar(theParName, convert(thePar), intE, theComment);
}

bool CXMLSwitchFileReader::writeBoolPar( const char* theParName,
	   bool thePar, const char *theComment)
{
	return writePar(theParName, convertbool(thePar), boolE, theComment);
}

bool CXMLSwitchFileReader::writeSwitchPar( const char* theParName,
	   int thePar, const char *theComment)
{
	char value[16];
	sprintf(value,"%d",thePar);
	return writePar(theParName, value, switchE, theComment);
}

bool CXMLSwitchFileReader::readStringPar( const char* theParName,
				   char* thePar, int maxLength)
{
	const XMLCh* val = readPar(theParName, stringE, thePar); 
	if(val==NULL)
		return false;

	DOMStringWrapper value = val;

	bool result = true;
	if(maxLength<(int)strlen(value))
	{
		strncpy(thePar,value,maxLength-1);
		thePar[maxLength-1] = '\0';
		cerr << "too long string parameter \"" << theParName << "\" - truncating\n";
	}
	else
		strcpy(thePar,value);
	return result;
}

bool CXMLSwitchFileReader::writePar( const char* theParName,
				   const char* thePar, const char *theComment)
{
	return writePar(theParName, thePar, stringE, theComment);
}


