#pragma once

#include "FVSystemOutputBase.h"

class FVVectorOutput : public FVSystemOutputBase
{
public:
	// Base
    static InputParameters validParams();
	FVVectorOutput(const InputParameters & parameters);
	
	// Main Output
	void output();

protected:
    FVVectorTupleMap extractVariableVector(std::string variable_name, std::string component_name, NumericVector<Real> & vector);
    void outputVector(std::string variable_name, std::string vector_name, std::string component_name, NumericVector<Real> & vector);

    std::string _vector_type;
	
};