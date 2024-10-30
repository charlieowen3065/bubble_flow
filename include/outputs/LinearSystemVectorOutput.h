#pragma once
#include "LinearSystemBaseOutput.h"

class LinearSystemVectorOutput : public LinearSystemBaseOutput
{
public:
    static InputParameters validParams();
    LinearSystemVectorOutput(const InputParameters & parameters);

    void output();
protected:
    void loadVector();
    void outputVectorOfSingleVariable(std::string variable_name);
    Real** createVectorOutputArrayOfSingleVariable(dof_id_type** component_array, int number_components);
    void variableVectors();
    void fullVector();

    NumericVector<Real> & _system_vector;
    std::string _vector_type; // Solution or Resiudal
};