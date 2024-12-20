#include "FVVectorOutput.h"

registerMooseObjectAliased("BubbleFlowApp", FVVectorOutput, "FVVectorOutput");

InputParameters
FVVectorOutput::validParams()
{
	InputParameters params = FileOutput::validParams();
}

FVVectorOutput::FVVectorOutput(const InputParameters & parameters)
  : FileOutput(parameters)
{

}

void
FVVectorOutput::output()
{

}