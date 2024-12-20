#pragma once

#include "FileOutput.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "FEProblemBase.h"

#include <chrono>
#include <ctime>

class FVVectorOutput : public FileOutput
{
public:
    static InputParameters validParams();
	FVVectorOutput(const InputParameters & parameters);
	
	// Main Output
	void output();

protected:

};