#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#include <assert.h>
#include <stdio.h>
#include <execinfo.h>
#include <unistd.h>

#define BACKTRACE_SIZE 10



inline void print_stack_trace() {
    void *array[10];
    size_t size;
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    // print out all the frames to stderr
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}

#ifdef _DEBUG

extern void debugBreakpoint();

#define ASSERT(f) {if(!(f))\
{                          \
    print_stack_trace();   \
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
