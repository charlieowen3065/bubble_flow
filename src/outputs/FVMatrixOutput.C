#include "FVMatrixOutput.h"

registerMooseObjectAliased("BubbleFlowApp", FVMatrixOutput, "FVMatrixOutput");

InputParameters
FVMatrixOutput::validParams()
{
    InputParameters params = FVSystemOutputBase::validParams();

    params.addParam<std::string>("matrix_type", "full", "Determines if the full jacobain is saved, or a variable-based jacobain {full, semi}");
    params.addParam<std::string>("output_format", "txt", "Determines the output format of the jacobain {txt, csv}");

    return params;
}

FVMatrixOutput::FVMatrixOutput(const InputParameters & parameters)
  : FVSystemOutputBase(parameters),
    _matrix_type(getParam<std::string>("matrix_type")),
    _output_format(getParam<std::string>("output_format"))
{
}

void
FVMatrixOutput::output()
{
    // Setup
    setupMatrix();

    // Output the matrix
    if (_matrix_type == "full")
        outputFullMatrix();
    else if (_matrix_type == "semi")
        outputSemiMatrix();
    else
        mooseError("Error: please enter valid jacobain output type ('" + _matrix_type + "' is invalid)");
}

void
FVMatrixOutput::setupMatrix()
{
    // Assemble matrix
	// Petsc stuff
	auto & petsc_options = _problem_ptr->getPetscOptions();
	auto & pars = _problem_ptr->solverParams();
	Moose::PetscSupport::petscSetOptions(petsc_options, pars);
	ExecFlagType flag = _current_execute_flag;
	_problem_ptr->execute(flag);
}

void
FVMatrixOutput::outputFullMatrix()
{
    auto & jacobian = static_cast<ImplicitSystem &>(_nl.system()).get_system_matrix();

    if (_output_format == "txt")
        saveJacobianToTextFile(jacobian, "full");
    else if (_output_format == "csv")
        saveJacobianToCSVFile(jacobian, "full");
    else
        mooseError("Error: please enter valid output format type ('" + _output_format + "' is invalid)");
}

void
FVMatrixOutput::outputSemiMatrix()
{
    // DOF Setup
    setDofMaps();

    // Collects the variable names
    std::vector<VariableName> variable_names;
    if ((_output_variable_name != "all") && (_output_variable_name != "none"))
        variable_names.push_back(_output_variable_name);
    else if (_output_variable_name == "all") 
        for (auto v : _nl.getVariableNames())
            variable_names.push_back(v);
    

    // TODO
}

void
FVMatrixOutput::saveJacobianToTextFile(SparseMatrix<Real> & jacobian, std::string variable_name)
{
	std::string base_name = "jacobian_" + variable_name;
	std::string filename = getFilename(base_name) + ".txt";
	
	// Create file
	std::ofstream txt_file;
	txt_file.open(filename);
	
	//std::cout << "HERE2" << std::endl;
	
	// Save Jacobian to the txt file
	jacobian.print(txt_file, true);
	
	// Close the txt file
	txt_file.close();
} 

void
FVMatrixOutput::saveJacobianToCSVFile(SparseMatrix<Real> & jacobian, std::string variable_name)
{
	std::string base_name = "jacobian_" + variable_name;
	std::string filename = getFilename(base_name) + ".csv";
	
	// Create file
	std::ofstream csv_file;
	csv_file.open(filename);

    // TODO
	
}