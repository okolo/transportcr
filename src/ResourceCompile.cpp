// ResourceComile.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include "base64.h"

using namespace std;

void printUsage(const char* progName)
{
	cerr << progName << "<identifier> <input file>\n output is written to stdout" << endl;
}

//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		printUsage(argv[0]);
		return 1;
	}
	const char* identifier = argv[1];

	fstream in;
	try{
		in.open(argv[2], ios_base::in);
	}catch(...)
	{
		cerr << "Failed to open file \"" << argv[2] << "\"" << endl;
		return 1;
	}
	unsigned int maxSize = 1024;
	char* buffer = new char[maxSize];
	maxSize--;
	unsigned int l = 0;
	int nLines = 0;
	for(in.read(buffer, maxSize),l=in.gcount(); l>0 ;in.read(buffer, maxSize),l=in.gcount())
	{
		if(nLines)
			cout << ",";
		else
			cout << "const std::string " << identifier << "_array[] = {\n";
		nLines++;
		buffer[l] = '\0';
		string encoded = base64_encode((unsigned char*)buffer, l);
		cout << "\"" << encoded << "\"\n";
		//test
		//std::string decoded = base64_decode(encoded);
		//cout << endl << endl << decoded << endl;
	}
	if(nLines>0)
	{
		cout << "};\nStringResource " << identifier << "(" 
			<< identifier << "_array, " << nLines << ");" << endl;
	}
	else
	{
		cout << "StringResource " << identifier << "(NULL,0);" << endl;	
	}
	delete[] buffer;
	return 0;
}

