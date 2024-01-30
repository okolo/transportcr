#include "Log.h"
#include "Addfunc.h"
#include <iostream>

Log logger;
const unsigned int Log::fMaxNumberOfIdenticalMessages = 10;

Log::Log(int aMaxNumberOfRecords, enum LogDetailLevel aDetailLevel):
fMaxNumberOfRecords(aMaxNumberOfRecords),
//fDetailLevel(aDetailLevel),
fCurNumberOfRecords(0)
{
}

void Log::Write(std::string aMessage, std::string aFile, int aLine, const LogDetailLevel aDetailLevel, bool aOnce)
{
	if(fCurNumberOfRecords > fMaxNumberOfRecords)
		return;
	std::ostringstream logStr;
	logStr << aFile << '(' << aLine << ")\t" << aMessage;
	Write(logStr.str(), aDetailLevel, aOnce);
}

void Log::Write(std::string aMessage, const LogDetailLevel aDetailLevel, bool aOnce)
{
	if(fCurNumberOfRecords > fMaxNumberOfRecords)
		return;
	const char* headers = "IMWE";//Info,Message,Warning,Error: the order is important!
	std::string msg;
	msg += headers[aDetailLevel];
	msg += ": ";
	msg += aMessage;
	///limit number of identical messages
	int hash = hash_string(msg.c_str());
	TInt& no = fSentMessages[hash];
	no.value++;
	if(no.value > fMaxNumberOfIdenticalMessages || (aOnce && no.value>1))
		return;
	if(aOnce)
		no.value = fMaxNumberOfIdenticalMessages;
	if(no.value == fMaxNumberOfIdenticalMessages)
		msg = "!!!LAST!!! " + msg;
	cout << msg << endl;
	fCurNumberOfRecords++;
	if(fCurNumberOfRecords > fMaxNumberOfRecords)
		cout << "!!! Maximal number of records achieved" << endl;
}

Log::~Log(void)
{
}
