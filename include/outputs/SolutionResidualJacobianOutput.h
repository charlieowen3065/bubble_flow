#pragma once

#include "FileOutput.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "FEProblemBase.h"

#include <chrono>
#include <ctime>

class SolutionResidualJacobianOutput : public FileOutput
{
public:
	static InputParameters validParams();
	SolutionResidualJacobianOutput(const InputParameters & parameters);
	
	// Main Output
	void output();
	
	// ---- Misc. Output Methods ---- //
	// CSV-Output
	void outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name);
	void outputToCSV(std::vector<std::vector<dof_id_type>> arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name);
	std::string getFilename(std::string BASE);  // Filename
	
	// DOF and Node getters
	dof_id_type** getDofIndices(std::vector<dof_id_type> & di, const unsigned int sys_number, const unsigned int var_number);
	void getDOFMapping(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping);
	std::vector<Real> get_x_y_z_id(const Point & p, const Real & id);
	void outputDOFMap();
	
	void displayTime();  // Time display
	Real computeL2Norm(std::vector<Real> residual_vector);  // L2-Norm
	
	// ---- Solution and Residual Methods ---- //
	Real** createSolutionResidualOutputArray(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping);
	Real** createDeltaSolutionOutputArray(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping);
	void outputSolutionResidual();
	void outputDeltaSolution();
	
	std::vector<Real> getResidualVector_SingleComponent(int component_number, std::vector<std::vector<dof_id_type>> & dof_mapping);
	std::vector<Real> getSolutionVector_SingleComponent(int component_number, std::vector<std::vector<dof_id_type>> & dof_mapping);

	void computeLinearIterationOutputs();
	void allLinearIterations(int number_linear_iterations);
	void getPreviousVectors();

	// ---- Jacobian Methods ---- //
	void JacobianOutput();
	void saveJacobianToTextFile(SparseMatrix<Real> & jacobian);
	void saveJacobianToCSVFile(SparseMatrix<Real> & jacobian);
	
	// ---- Debug Methods ---- //
	void debugOutput();
	
protected:
	
	// ---------- System Variables ---------- //
	
	// Nonlinear System Variables
	THREAD_ID _tid;
	const unsigned int _nl_sys_num;  // The nonlinear system number we should output degree of freedom information for
	NonlinearSystemBase & _nl;
	
	// System Variables
	const System & _sys;
	const DofMap & _dof_map;
	const unsigned int _sys_number;
	std::string _system_name;  // The name of the system to extract DOF information

	// Meshing
	MooseMesh & _mesh;  // Reference to the mesh object
	dof_id_type _num_nodes;
	
	// Nonlinear-Variable Information
	std::string _output_variable_name;
	MooseVariableFieldBase & _output_variable;
	int _output_variable_number;
	unsigned int _number_componets;
	bool _is_vector;
	
	// Solution, Residual, and Jacobian References
	NumericVector<Real> & _solution;
	NumericVector<Real> & _residual;
	NumericVector<Real> & residual_full;
	
	std::vector<Real> _prev_it_solution;
	std::vector<Real> _prev_it_residual;
	
	
	int _linear_itt_counter;
	bool _copy_var;
	
	// ---------- User-Inputs ---------- //
	
	// Solution (& DeltaS), Residual, Jacobian Inputs
	std::vector<std::string> _which_outputs;
	bool _output_solution;
	bool _output_delta_solution;
	bool _output_residual;
	bool _output_jacobian;
	bool _output_DOF_map;
	
	// Timing Inputs
	std::string _which_time;
	bool _show_converged_time;
	bool _show_nl_time;
	bool _show_l_time;
	std::string _time_units;
	int START_TIME_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	
	// Solution & Residual Specific Inputs
	bool _solution_residual_seperate_files;
	
	// Jacobian-Specific Inputs
	bool _display_jacobian;
	bool _save_jacobian_txt;
	bool _save_jacobian_csv;
	
	// Debug
	bool _SRJ_debug;
	bool _show_L2_norm;
	
	
	
};