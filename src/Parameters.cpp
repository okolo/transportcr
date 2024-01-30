#include "Parameters.h"

IParReader* Parameters::fReader = 0;

IParWriter* IParWriter::fInstance = 0;

Parameters::Parameters()
{
	if(!fReader)
		ThrowError("Attempt to create parameter class before initializing Parameter Reader");

}

void Parameters::SetReader(IParReader* aReader)
{
	fReader = aReader;
}
