#pragma once
#include "LinearSystemBaseOutput.h"

class SolutionVectorOutput : public LinearSystemBaseOutput
{
public:
    static InputParameters validParams();
    SolutionVectorOutput(const InputParameters & parameters);

    void output();
protected:
    void outputSolutionOfSingleVariable(std::string variable_name);
    Real** createSolutionOutputArrayOfSingleVariable(dof_id_type** component_array, int number_components);
    void variableSolutions();
    void fullSolution();
};