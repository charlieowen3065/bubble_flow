#pragma once

#include "FVSystemOutputBase.h"

class FVMatrixOutput : public FVSystemOutputBase
{
public:
	// Base
    static InputParameters validParams();
	FVMatrixOutput(const InputParameters & parameters);
	
	// Main Output
	void output();

protected:
    // Matrix setup
    void setupMatrix();

    void outputFullMatrix();
    void outputSemiMatrix();

    void saveJacobianToTextFile(SparseMatrix<Real> & jacobian, std::string variable_name);
    void saveJacobianToCSVFile(SparseMatrix<Real> & jacobian, std::string variable_name);

    std::string _matrix_type;
    std::string _output_format;
	
};