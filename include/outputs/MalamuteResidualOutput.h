#pragma once

// MOOSE includes
#include "FileOutput.h"
#include "libmesh/dof_map.h"
#include "NonlinearSystem.h"

class MooseMesh;

class MalamuteResidualOutput : public FileOutput
{
public:
	static InputParameters validParams();

	MalamuteResidualOutput(const InputParameters & parameters);
	
	void output();
	
	void getDOFMapping(std::vector<std::string> & resdiual_names, std::vector<std::vector<dof_id_type>> & dof_mapping);
	dof_id_type** getDofIndices(std::vector<dof_id_type> & di, const unsigned int sys_number, const unsigned int var_number);
	
	Real** getResiduals(std::vector<std::string> & resdiual_names, std::vector<std::vector<dof_id_type>> & dof_mapping);
	
	void outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols);
	
	std::vector<Real> get_x_y_z_id(const Point & p, const Real & id);
	std::string getFilename();
	
protected:
	
	// System Parameters
	std::string _system_name;  // The name of the system to extract DOF information
	
	MooseMesh & _mesh;  // Reference to the mesh object
	dof_id_type _num_nodes;
	
	THREAD_ID _tid;
	
	const unsigned int _nl_sys_num;  // The nonlinear system number we should output degree of freedom information for
	NonlinearSystemBase & _nl;
	NumericVector<Real> & _all_residuals;
	
	const System & _sys;
	const DofMap & _dof_map;
	const unsigned int _sys_number;
	
	// User Inputs
	std::string _output_variable;
	
	bool _write_file;  // Flag for controlling outputting console information to a file
	bool _write_screen;  // Flag for controlling outputting console information to screen
	
	bool _get_dof_map;
	bool _get_res_map;
	bool _check_has_var;
	MooseVariableFieldBase & _mal_var;
	
	int _mal_var_number;
	unsigned int _num_comps;
	bool _is_vector;
	

};