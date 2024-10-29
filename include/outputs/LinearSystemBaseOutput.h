#pragma once

#include "FileOutput.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "FEProblemBase.h"

#include <chrono>
#include <ctime>

class LinearSystemBaseOutput : public FileOutput
{
public:
    static InputParameters validParams();
    LinearSystemBaseOutput(const InputParameters & parameters);

    // Main Output
    void output();

protected:
    // Setup functions
    void setupTimeVariables();
    std::string getFilename(std::string BASE);
    void outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name);
    std::vector<Real> getNodeCoordinates(const Point & p, const Real & id);
    dof_id_type** getDOFIndices_SingleVariable(std::vector<dof_id_type> & di, const unsigned int sys_number,
                                               const unsigned int var_number, unsigned int number_componets);

    // -------- Class Variables -------- //
    /// Class Initilization
    // Nonlinear system variables
    THREAD_ID _tid;
    const unsigned int _nl_sys_num;
	NonlinearSystemBase & _nl;
    
    // Nonlinear-Variable Information
	std::string _output_variable_name;
    
    // System variables
    const System & _sys;
	const DofMap & _dof_map;
	const unsigned int _sys_number;
	std::string _system_name;
    
    // Meshing
    MooseMesh & _mesh;
	dof_id_type _num_nodes;
    
    // Timing Inputs
    std::string _which_time;
	std::string _time_units;

    /// Latter Initilizations
    // setupTimeVariables
    bool _show_converged_time;
	bool _show_nl_time;
	bool _show_l_time;

    // START-TIME
    int START_TIME_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
};