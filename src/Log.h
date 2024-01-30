#ifndef Log_H_ALREADY_INCLUDED
#define Log_H_ALREADY_INCLUDED

#include <string>
#include <map>

class Log
{
public:
	enum LogDetailLevel
	{//the order is important!
		EVerbose = 0,//use for debug info output
		EMessage,
		EWarning,
		EError
	};
	Log(int aMaxNumberOfRecords = 100000, const LogDetailLevel aDetailLevel = EError);
	~Log(void);
	void Write(std::string aMessage, std::string aFile, int aLine, const LogDetailLevel aDetailLevel, bool aOnce=false);
	void Write(std::string aMessage, const LogDetailLevel aDetailLevel = EMessage, bool aOnce=false);
private:
	class TInt
	{
	public:
		TInt() : value(0){}
		unsigned int value;
	};
	int							fMaxNumberOfRecords;
	//const						LogDetailLevel fDetailLevel;
	int							fCurNumberOfRecords;
	static const unsigned int	fMaxNumberOfIdenticalMessages;
	std::map<int, TInt>			fSentMessages;
};

extern Log logger;

#define LOG_MESSAGE(par1) { std::ostringstream __logStr;\
	__logStr << (par1); \
	logger.Write(__logStr.str(), Log::EMessage); }

#define LOG_MESSAGE2(par1, par2) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2); \
	logger.Write(__logStr.str(), Log::EMessage); }

#define LOG_MESSAGE3(par1, par2, par3) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3); \
	logger.Write(__logStr.str(), Log::EMessage); }

#define LOG_MESSAGE4(par1, par2, par3, par4) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3) << (par4); \
	logger.Write(__logStr.str(), Log::EMessage); }

#define LOG_WARNING(par1) { std::ostringstream __logStr;\
	__logStr << (par1); \
	logger.Write(__logStr.str(), Log::EWarning); }

#define LOG_WARNING2(par1, par2) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2); \
	logger.Write(__logStr.str(), Log::EWarning); }

#define LOG_WARNING3(par1, par2, par3) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3); \
	logger.Write(__logStr.str(), Log::EWarning); }

#define LOG_WARNING4(par1, par2, par3, par4) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3) << (par4); \
	logger.Write(__logStr.str(), Log::EWarning); }

#define LOG_ERROR(par1) { std::ostringstream __logStr;\
	__logStr << (par1); \
	logger.Write(__logStr.str(), Log::EError); }

#define LOG_ERROR2(par1, par2) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2); \
	logger.Write(__logStr.str(), Log::EError); }

#define LOG_ERROR3(par1, par2, par3) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3); \
	logger.Write(__logStr.str(), Log::EError); }

#define LOG_ERROR4(par1, par2, par3, par4) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3) << (par4); \
	logger.Write(__logStr.str(), Log::EError); }

#define LOG_ERROR_ONCE(par1) { std::ostringstream __logStr;\
	__logStr << (par1); \
	logger.Write(__logStr.str(), Log::EError, true); }

#define LOG_ERROR2_ONCE(par1, par2) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2); \
	logger.Write(__logStr.str(), Log::EError, true); }

#define LOG_ERROR3_ONCE(par1, par2, par3) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3); \
	logger.Write(__logStr.str(), Log::EError, true); }

#define LOG_ERROR4_ONCE(par1, par2, par3, par4) { std::ostringstream __logStr;\
	__logStr << (par1) << (par2) << (par3) << (par4); \
	logger.Write(__logStr.str(), Log::EError, true); }

#endif
