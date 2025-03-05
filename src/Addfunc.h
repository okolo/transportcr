#if !defined(ADDFUNC_H_INCLUDED)
#define ADDFUNC_H_INCLUDED


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>

#ifdef USE_GSL
#include <gsl/gsl_vector_double.h>
#endif

#include "debug.h"
#include "Log.h"

using namespace std;

#define DIR_DELIMITER_CH '/'
#define DIR_DELIMITER_STR "/"


#define DATA_DIR "." DIR_DELIMITER_STR "tables" DIR_DELIMITER_STR
#define ALL_PLT_DIR "." DIR_DELIMITER_STR "results" DIR_DELIMITER_STR
#define PLT_DIR_LINK ALL_PLT_DIR "local_tmp"


#define STRING2(x) #x
#define STRING(x) STRING2(x)

extern const double Pi;

#ifndef INFINITY
#define INFINITY 1e308
#endif

#define UINT unsigned int
typedef double* PDouble;
typedef double** matrix;
typedef std::vector<double>	StdVector;

extern string plt_local_dir;
extern const char* plt_local_c;

enum EReturnCode
{
	retOk = 0,
	retUsageError = -1,
	retIOError = 10,
	//retWrongSettings = 11,
	retInvalidSettingsFile = 12,
	//retNotEnoughAccuracy = 13,
	//retCalcError = 14,
	//retNotEnoughMemory = 15,
	retBufferOverflow = 16,
	//retGeneralExeption = 17,
	//retTableInitError = 18,
	//retTestNormalExit = 19,
	//retWrongFuncParam = 20,
	//retAssertionFailed = 21,
	//retInvalidValue = 22
};

#define BACKTRACE_MAX_SIZE 20
extern void *backtrace_array[];
extern char **stack_trace_symbols;
extern size_t stack_trace_size;
extern void print_saved_stack_trace();
inline void save_stack_trace() {
    if (stack_trace_symbols)
        free(stack_trace_symbols); // clear saved info
    // get void*'s for all entries on the stack
    stack_trace_size = backtrace(backtrace_array, BACKTRACE_MAX_SIZE);
    stack_trace_symbols = backtrace_symbols(backtrace_array, stack_trace_size);
}
void ThrowError(std::string aMsg);//safely throw error (aMsg can be stack variable)
void MemoryExit(const char* str=NULL);
void OpenExit(string str);
void ReportError(const char* errStr);
FILE *Fopen(string _file,const char *_filetype);
FILE *FopenL(string _file,const char *_filetype);//throws _file
FILE *FopenL(string _file, string _dir, const char *_filetype);//throws _file
FILE *Fopen(string _file, string _dir, const char *_filetype, bool exitOnError=true);
int IsEqual(double,double);
void* Calloc(size_t count, size_t eltsize, const char* erStr=NULL);
void Exit(int code);
bool fileExists(std::string fileName);
bool directory_exists(std::string aPath);
unsigned int hash_string(const char * s);

int Mkdir(string _dir);

int CopyFileOrDir(std::string from, std::string to);//return nonzero if seccess
int MoveFileOrDir(std::string from, std::string  to);//return nonzero if seccess
int RemoveFileOrDir(std::string name);//return nonzero if seccess

const char* GetFileName(const char* path);


inline int Round(double a)
{
	return (int)(a+0.5);
}

void trim_left(std::string &str, const char* whiteSpaceChars = " \t\r\n");

void trim_right(std::string &str, const char* whiteSpaceChars = " \t\r\n");

inline void trim(std::string &str, const char* whiteSpaceChars = " \t\r\n")
{
	trim_left(str, whiteSpaceChars);
	trim_right(str, whiteSpaceChars);
}

template<class T> std::string ToString (const T& aValue)
{
	std::ostringstream logStr;
	logStr << aValue;
	return logStr.str();
}

///read simple types in binary and text modes
template<class T> bool read(ifstream& aStream, T& aValRef, bool aBinaryMode = false)
{
	if (aStream.eof()) {
		return false;
	}

	if (aBinaryMode) {
		T value = 0;
		aStream.read((char *) &value, sizeof(T));
		aValRef = value;
	}else{
		aStream >> aValRef;
	}

	if (aStream.fail()) {
		return false;
	}
	return true;
}

template<class T> void readL(ifstream& aStream, T& aValRef, bool aBinaryMode = false)
{
	if(!read(aStream, aValRef, aBinaryMode))
		throw aStream.eof() ? "unexpected end of file" : "data format error";
}

template<class T> bool write(ofstream& aStream, T aVal, unsigned char aDelimiter, bool aBinaryMode = false)
{
	if (aBinaryMode) {
		aStream.write((const char*)&aVal, sizeof(T));
	}else{
		aStream << aDelimiter << aVal;
	}
	return !aStream.fail();
}

template<class T> void writeL(ofstream& aStream, T aVal, unsigned char aDelimiter, bool aBinaryMode = false)
{
	if(!write(aStream, aVal, aDelimiter, aBinaryMode))
		throw "writing error";
}

template<class T> inline T MAX (const T& aVal1, const T& aVal2)
{
	return aVal1>aVal2 ? aVal1 : aVal2;
}

template<class T> inline T MIN (const T& aVal1, const T& aVal2)
{
	return aVal1<aVal2 ? aVal1 : aVal2;
}


//  additional functions used for debugging
int IsInfinity(double _val);

int IsValid(double _val);/* returns zero if _val<0 or is Infinity */
int IsValid(double _val, int _line, const char* _file, bool aThrowIfNaN = true); /* returns zero if _val<0 or is Infinity
writes to output information about error position in source file */

#ifdef _DEBUG
#define ASSERT_VALID_NO(f) ASSERT(IsValid(f))
#define RERURN_VALID(f) {double _ret=(f); ASSERT(IsValid(_ret)); return _ret;}
#define VERIFY_VALID_NO(f) ASSERT(IsValid(f))
#else
#define ASSERT_VALID_NO(f)          ((void)0)
#define RERURN_VALID(f)				return (f)
#define VERIFY_VALID_NO(f)          ((void)(f))
#endif

double Exp(double _val);

#define TESTVALUE(_val,_val1) if(!IsValid(_val,__LINE__,__FILE__))\
{\
   _val=_val1;\
}
#define WARN(warning) logger.Write(warning, __FILE__, __LINE__, Log::EWarning)
#define ERR(warning) logger.Write(warning, __FILE__, __LINE__, Log::EError)
#define WARN_ONCE(warning) logger.Write(warning, __FILE__, __LINE__, Log::EWarning, true)
#define ERR_ONCE(warning) logger.Write(warning, __FILE__, __LINE__, Log::EError, true)


template<class type>
class SafePtr
{
public:
	SafePtr(type* _pType=NULL):pType(_pType){};
	virtual ~SafePtr(){delete pType;};
	SafePtr& operator=(type* _pType){delete pType; pType=_pType; return (*this);};
	bool isNull() const {return (pType==NULL);};
	inline type* Detach() { type* result = pType; pType = 0; return result; }

	operator type*(){ return pType; };
	operator const type*() const { return pType; };

	//operator type&(){ASSERT(pType);return *pType;};
	//operator const type&() const {ASSERT(pType);return *pType;};

	type& operator*() {ASSERT(pType);return *pType;};
	const type& operator*() const {ASSERT(pType);return *pType;};

	inline type* operator->() const
	{
		ASSERT(pType);
		return pType;
	}
private:
	type* pType;
};

inline double average(double _part,double _val_1,double _val_2)
{
	return (1.0-_part)*_val_1+_part*_val_2;
} 

class TReferencedObj{
public:
	TReferencedObj():iRefCount(0){};
	void addRef(){
		iRefCount++;
	}
	void releaseRef(){
		iRefCount--;
		if (iRefCount<=0) {
			delete this;
		}
	}
	inline int RefCount() const{return iRefCount;};
	virtual ~TReferencedObj(){};
private:
	int iRefCount;
};

template <class T>//class T should have addRef() & releaseRef() methods
class SmartPtr
{
public:
  SmartPtr(T* pointee = NULL) : iPointee(pointee){
		if (iPointee) iPointee->addRef();
  };

  SmartPtr(const SmartPtr<T>& other) : iPointee(other.iPointee){
		if (iPointee) iPointee->addRef();
  };

  inline SmartPtr& operator=(T* pointee){
	  if(iPointee==pointee)
		  return *this;
	  if (iPointee) iPointee->releaseRef();
	  iPointee = pointee;
	  if (pointee) iPointee->addRef();
	  return *this;
  }

  inline SmartPtr& operator=(const SmartPtr<T>& other){
	  if (iPointee) iPointee->releaseRef();
	  iPointee = other.iPointee;
	  if (iPointee) iPointee->addRef();
	  return *this;
  }

  ~SmartPtr(){
		if (iPointee) iPointee->releaseRef();
  }

  inline bool operator==(T* pointee) const{
	  return iPointee==pointee;
  }

  inline T& operator*() const
  {
	return *iPointee;
  }

  inline T* operator->() const
  {
	return iPointee;
  }
  
  inline operator T*() const
  {
	return iPointee;
  }

  inline T& operator[](int) const{
	  //! although operator T* is defined, smart pointers should not be used as C-arrays
	  NOT_IMPLEMENTED;
	  return *iPointee;
  }

private:
  T* iPointee;
};



#endif// end of file
