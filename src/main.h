#if !defined(E_PROPAG_H_INCLUDED)
#define E_PROPAG_H_INCLUDED

#include "Parameters.h"

class Application : Parameters
{
public:
	static int main(int argc,char *argv[]);
private:
	Application();
	void run();
	class CInjectionSpectra* CreateSource();

	static void GslFrrorHandler (const char * reason,
	        const char * file,
	        int line,
	        int gsl_errno);
	static void usageExit();
	static int UnitTest(int argc,char *argv[]);

	static const char* progName;
};


#endif //end of file
