#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#include <assert.h>

#ifdef _DEBUG

extern void debugBreakpoint();

#define ASSERT(f) {if(!(f))\
{\
	debugBreakpoint();\
	assert(f);\
};}// assert(f) will print information about condition <f> and place of error
/// debugBreakpoint method may be used in debugging (breakpoint must be set
/// in this method)


#define VERIFY(f) ASSERT(f)
#define DEBUG_ONLY(f)      (f)

#else   // _DEBUG

#define ASSERT(f)          ((void)0)
#define VERIFY(f)          ((void)(f))
#define DEBUG_ONLY(f)      ((void)0)

#endif  // _DEBUG

#define NOT_IMPLEMENTED	ASSERT(false);//not implemented yet

#endif//end of file
