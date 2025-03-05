
#include <stdlib.h>
#include <math.h>
#include "Addfunc.h"
#include <string.h>
#include <fstream>
#include <dirent.h>
#include <stdio.h>

void *backtrace_array[BACKTRACE_MAX_SIZE];
char **stack_trace_symbols = NULL;
size_t stack_trace_size = 0;

void print_saved_stack_trace() {
    fprintf(stderr, "Stack trace:\n");
    for (size_t i = 0; i < stack_trace_size; i++) {
        fprintf(stderr, "%s\n", stack_trace_symbols[i]);  // Function names may not appear without -rdynamic
    }
    free(stack_trace_symbols);
    stack_trace_symbols = NULL;
    stack_trace_size = 0;
}

using namespace std;

const double Pi=M_PI;
static std::string LastErrorMessage;

#include <algorithm>

#ifdef USE_GSL
template<> SafePtr<gsl_vector>::~SafePtr()
{
	if(pType)
		gsl_vector_free(pType);
}
#endif

static std::string errorBuffer;

void ThrowError(std::string aMsg)
{
    errorBuffer = aMsg;
    save_stack_trace();
	throw errorBuffer.c_str();
}

const char* LastError()
{
	return errorBuffer.c_str();
}

void trim_left(std::string &str, const char* whiteSpaceChars)
{
	for(int i=0; whiteSpaceChars[i]!='\0'; i++)
		str.erase(0, str.find_first_not_of(whiteSpaceChars[i]));
}

void trim_right(std::string &str, const char* whiteSpaceChars)
{
	for(int i=0; whiteSpaceChars[i]!='\0'; i++)
		str.erase(str.find_last_not_of(whiteSpaceChars[i]) + 1, std::string::npos);
}

unsigned int hash_string(const char * s)
{
    unsigned int hash = 0;

    for(; *s; ++s)
    {
        hash += *s;
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}

bool fileExists(std::string fileName)
{
	ifstream f(fileName.c_str());
	bool result = (bool)f.is_open();
	if(result)
		f.close();
	return result;
}

bool directory_exists(std::string aPath)
{
	DIR *pdir=opendir (aPath.c_str());
	bool result = (pdir!=0);
	if(result)
		closedir(pdir);
	return result;
}

void MemoryExit(const char* str)
{
	if(str!=NULL)
		ThrowError(ToString("not enough memory ") + str);
	else
		ThrowError("not enough memory");
};

void* Calloc(size_t count, size_t eltsize,const char* errStr)
{
	void* ptr=calloc(count,eltsize);
	if(ptr==NULL)
		MemoryExit(errStr);
	return ptr;
}

void ReportError(const char* errStr)
{
	cerr << errStr << endl;
#ifdef _DEBUG
	debugBreakpoint();
#endif
}

void OpenExit(string str)
{
	LastErrorMessage = "can't open file " + str;
	LOG_ERROR(LastErrorMessage);
	Exit(retIOError);
};

FILE *Fopen(string _file,const char *_filetype)
{
	FILE *out=fopen(_file.c_str(), _filetype);
	if(out==NULL)
		OpenExit(_file);
	return out;
};

FILE *FopenL(string _file, const char *_filetype)
{
	FILE *out=fopen(_file.c_str(),_filetype);
	if(out==NULL){
		ThrowError("can't open file " + _file);
	}
	return out;
};

FILE *FopenL(string _file, string _dir, const char *_filetype)//throws _file
{
	string path(_dir);
	if(_dir[_dir.length()-1]!=DIR_DELIMITER_CH)
		path += DIR_DELIMITER_CH;

	path += _file;

	FILE *out=fopen(path.c_str(),_filetype);
	if(out==NULL){
		ThrowError("failed to open '" + path + "' as '" + _filetype + "'");
	}
	return out;
}

FILE *Fopen(string _file, string _dir, const char *_filetype, bool exitOnError)
{
	FILE *out = NULL;
	try{
		out = FopenL(_file, _dir, _filetype);
	}catch(const char*){
	if((out==NULL)&&exitOnError)
		Exit(retIOError);
	}
	return out;
}

int IsEqual(double _val1,double _val2)
{
	const double accuracy=1e-7;
	if(_val1==_val2)
		return 1;
	if ((_val1+_val2)==0)
		return 0;
	if(fabs((_val1-_val2)/(_val1+_val2))<accuracy)
		return 1;
	return 0;
}
//functions need for debugging

int IsInfinity(double _val) /* returns nonzero value if argument is Infinity */
{
	if(_val>1.6e308)
		return 1;
	if(_val<(-1.6e308))
		return -1;
	return 0;
}

int IsValid(double _val, int _line, const char* _file, bool aThrowIfNaN) /* returns zero if _val<0 or is Infinity */
{
	if(_val>=0 && _val < 1.6e308)
	{
		return 1;
	}
	string str;
	if(_val < 0)
		str = "negative value";
	else
	{
		str = "NaN value";
		if(aThrowIfNaN)
		{
			logger.Write(str, _file, _line, Log::EError);
			throw "Program aborted due to assertion failure";
		}
	}
#ifdef EXIT_ON_INVALID_VAL
	logger.Write(str, _file, _line, Log::EError);
	Exit(retInvalidValue);
#else
	logger.Write(str, _file, _line, Log::EWarning);
#endif
	return 0;
}

int IsValid(double _val) /* returns zero if _val<0 or is Infinity */
{
	if(_val>=0 && _val < 1.6e308)
	{
		return 1;
	}
	return 0;
}

double Exp(double _val)
{
	if(_val>709.666)
	{
		LOG_ERROR3("exp(", _val, ") > 1.6e308 floating point overflow")
	}
	if(_val<(-709.666))
		return 0.0;
	return exp(_val);
}

void Exit(int code)
{
	exit(code);
}

int Mkdir(string _dir)
{
	string command = "mkdir ";

	int length=_dir.length();
	char ch=_dir[length-1];
	if(ch==DIR_DELIMITER_CH)
	{
		string dir(_dir.c_str(),length-1);
		command = command + dir;
	}
	else
		command = command + _dir;
	return system(command.c_str());
}

/*void MkPltdir(const char *_dir)
{
	Mkdir((plt_local_dir + _dir).c_str());
}*/

int CopyFileOrDir(std::string from, std::string to)//return nonzero if seccess
{
	std::string command = "cp -R \"" + from + "\" \"" + to + "\"";
	return !system(command.c_str());
}

int MoveFileOrDir(std::string from, std::string to)//return nonzero if seccess
{
    std::string command = "mv -f \"" + from + "\" \"" + to + "\"";
    return !system(command.c_str());
}

int RemoveFileOrDir(std::string name)//return nonzero if seccess
{
	int result=!remove(name.c_str());
	if(!result)
	{
		result=!system(("rm -rf \"" + name + "\"").c_str());
	}
	return result;
}

const char* GetFileName(const char* path)
{
	const char* str1=strrchr(path,'/');//searching for last character "/" in the string
	const char* str2=strrchr(path,'\\');//searching for last character "\" in the string
	if(str2>str1)
		str1=str2;
	if(str1==NULL)
		str1=path;//input parameter contains only file name without full path
	else
		str1++;
	return str1;
}

#ifdef _DEBUG

void debugBreakpoint()
{
	printf("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	"Debug breakpoint has been reached\n"
	"Use debugger and set breakpoint on method debugBreakpoint()\n"
	"to catch error occured\n");
}

#endif

