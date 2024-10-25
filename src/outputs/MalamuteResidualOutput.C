// MOOSE includes
#include "MalamuteResidualOutput.h"
#include "FEProblem.h"
#include "KernelBase.h"
#include "MooseApp.h"
#include "Moose.h"
#include "Conversion.h"
#include "MooseMesh.h"


#include "libmesh/fe.h"

// compiler includes (for type demangling)
#include <cxxabi.h>
#include <fstream>

registerMooseObjectAliased("MalamuteApp", MalamuteResidualOutput, "MalamuteResidualOutput");

InputParameters
MalamuteResidualOutput::validParams()
{
  // Get the parameters from the base class
  InputParameters params = FileOutput::validParams();

  // Screen and file output toggles
  params.addParam<bool>("output_screen", false, "Output to the screen");
  params.addParam<bool>("output_file", true, "Output to the file");
  params.addParam<std::string>("system_name", "nl0", "System to output");

  // By default this only executes on the initial timestep
  params.set<ExecFlagEnum>("execute_on", true) = EXEC_INITIAL;

  params.addParam<NonlinearSystemName>(
      "nl_sys", "nl0", "The nonlinear system that we should output information for.");
  
  // New Methods
  params.addParam<bool>("get_dof_map", false, "Output the DOF Map");
  params.addParam<bool>("get_residuals_map", false, "Output the Residuals Map");
  params.addParam<bool>("check_has_var", "Checks to see if the problem has the variable");
  params.addParam<std::string>("output_variable", "Checks to see if the problem has a variable");
  
  
  return params;
}

MalamuteResidualOutput::MalamuteResidualOutput(const InputParameters & parameters)
  : FileOutput(parameters),
    
	// System Parameters
	_system_name(getParam<std::string>("system_name")),
	
	_mesh(_problem_ptr->mesh()),
	_num_nodes(_mesh.nNodes()),
	
	_tid(getParam<THREAD_ID>("_tid")),
	
    _nl_sys_num(_problem_ptr->nlSysNum(getParam<NonlinearSystemName>("nl_sys"))),
	_nl(_problem_ptr->getNonlinearSystemBase(_nl_sys_num)),
	_all_residuals(_nl.RHS()),
	
	_sys(_problem_ptr->es().get_system(getParam<std::string>("system_name"))),
	_dof_map(_sys.get_dof_map()),
	_sys_number(_dof_map.sys_number()),
	
	// User Inputs
	_output_variable(getParam<std::string>("output_variable")),
	
	_write_file(getParam<bool>("output_file")),
    _write_screen(getParam<bool>("output_screen")),
	
	_get_dof_map(getParam<bool>("get_dof_map")),
	_get_res_map(getParam<bool>("get_residuals_map")),
	_check_has_var(getParam<bool>("check_has_var")),
	_mal_var(_nl.getVariable(_tid, _output_variable))
	
	
{
	if (_check_has_var){
		if (_problem_ptr->hasVariable(_output_variable)){
			_console << _output_variable << " : True " << std::endl;
		} else {
			_console << _output_variable << " : False " << std::endl;
		}
	}
	
	
	
}

void
MalamuteResidualOutput::output()
{
	_mal_var_number = _mal_var.number();
	_num_comps = _mesh.nodePtr(0)->n_comp(_sys_number, _mal_var_number);
	if (_num_comps == 1) {_is_vector = false;} else {_is_vector = true;}
		
	std::vector<std::string> resdiual_names;
	std::vector<std::vector<dof_id_type>> dof_mapping;
	getDOFMapping(resdiual_names, dof_mapping);
	
	Real** residuals = getResiduals(resdiual_names, dof_mapping);
	
	std::vector<std::string> column_names;
	column_names.push_back("id");
	column_names.push_back("x");
	column_names.push_back("y");
	column_names.push_back("z");
	for (auto r_name : resdiual_names) {
		column_names.push_back(r_name);
	}
	
	outputToCSV(residuals, column_names, _num_nodes, 4 + _num_comps);

}

dof_id_type**
MalamuteResidualOutput::getDofIndices(std::vector<dof_id_type> & di, const unsigned int sys_number, const unsigned int var_number)
{
	
	// Declare return array
	dof_id_type** component_array = new dof_id_type*[_num_comps];
	
	
	// Get dof indices
	for (auto comp_number : make_range(_num_comps))
	{
		component_array[comp_number] = new dof_id_type[_num_nodes];
		
		for (auto node_idx : make_range(_num_nodes))
		{
			const Node * node_i = _mesh.nodePtr(node_idx);
			
			dof_id_type dof_i = node_i->dof_number(_sys_number, _mal_var_number, comp_number);
			
			di.push_back(dof_i);
			component_array[comp_number][node_idx] = dof_i;
			
		}
		
	}
	
	return component_array;
}

void
MalamuteResidualOutput::getDOFMapping(std::vector<std::string> & resdiual_names, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	
	std::vector<dof_id_type> dof_indices_full;
	dof_id_type** component_array = getDofIndices(dof_indices_full, _sys_number, _mal_var_number);
		
	// Print DOF Mapping
	if (_write_screen and _get_dof_map)
	{
		std::cout << "\n DOF MAPPING \n" << std::endl;
		std::vector<std::string> dof_mapping_str;
		if (_is_vector == false) {
			std::string dof_str = _output_variable + ": [";
			for (auto i : dof_indices_full) {
				dof_str += std::to_string(i) + " ";
			}
			dof_str += "]";
			
			dof_mapping_str.push_back(dof_str);	
		} else {
			
			std::vector<std::string> vec_comps = {"_x", "_y", "_z"};
			for (auto c_i : make_range(_num_comps))
			{
				std::string dof_str = _output_variable + vec_comps[c_i] + ": [";
				for (auto i : make_range(_num_nodes)) 
				{
					dof_str += std::to_string(component_array[c_i][i]) + " ";
				}
				dof_str += "]";
			
				dof_mapping_str.push_back(dof_str);	
			}	
		}
		
		for (auto top_i : dof_mapping_str){
			std::cout << top_i << "\n" << std::endl;
		}
	}
	
	// Place DOF Mapping into vector
	if (_is_vector == false) {
		resdiual_names.push_back(_output_variable);
		std::vector<dof_id_type> dof_map_i;
		for (auto i : dof_indices_full) {
			dof_map_i.push_back(i);
		}
		
		dof_mapping.push_back(dof_map_i);	
	} else {
		
		std::vector<std::string> vec_comps = {"_x", "_y", "_z"};
		for (auto c_i : make_range(_num_comps))
		{
			resdiual_names.push_back(_output_variable + vec_comps[c_i]);
			std::vector<dof_id_type> dof_map_i;
			for (auto i : make_range(_num_nodes)) 
			{
				dof_map_i.push_back(component_array[c_i][i]);
			}
		
			dof_mapping.push_back(dof_map_i);		
		}	
	}

}

Real**
MalamuteResidualOutput::getResiduals(std::vector<std::string> & resdiual_names, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	
		
	// Print Residuals Mapping
	if (_write_screen and _get_res_map)
	{
		std::cout << "\n RESIDUALS \n" << std::endl;
		std::vector<std::string> residual_strs;
		if (_is_vector == false) {
			std::stringstream res_str;
			res_str<<resdiual_names[0]<<": [";
			for (auto i : make_range(_num_nodes)) {
				dof_id_type dof_i = dof_mapping[0][i];
				Real res_i = _all_residuals(dof_i);
				res_str<<res_i<<" ";
			}
			res_str<<"]";
			
			residual_strs.push_back(res_str.str());	
		} else {
			for (auto c_i : make_range(_num_comps))
			{			
				std::stringstream res_str;
				res_str<<resdiual_names[c_i]<<": [";
				for (auto i : make_range(_num_nodes)) 
				{
					dof_id_type dof_i = dof_mapping[c_i][i];
					Real res_i = _all_residuals(dof_i);
					res_str<<res_i<<" ";
				}
				res_str<<"]";
				
				residual_strs.push_back(res_str.str());
			}	
		}
		
		
		for (auto top_i : residual_strs){
			std::cout << top_i << "\n" << std::endl;
		}
	}
	
	// Place Residuals into array
	Real** residuals = new Real*[_num_nodes];
	
	for (dof_id_type node_idx : make_range(_num_nodes))
	{
		const Node * node_i = _mesh.nodePtr(node_idx);
		std::vector<Real> x_y_z_id = get_x_y_z_id(*node_i, node_idx);
		
		residuals[node_idx] = new Real[4 + _num_comps];
		residuals[node_idx][0] = x_y_z_id[3];  // id
		residuals[node_idx][1] = x_y_z_id[0];  // x
		residuals[node_idx][2] = x_y_z_id[1];  // y
		residuals[node_idx][3] = x_y_z_id[2];  // z
		
		for (auto c_i : make_range(_num_comps))
		{
			dof_id_type dof_i = dof_mapping[c_i][node_idx];
			Real res_i = _all_residuals(dof_i);
			
			residuals[node_idx][4 + c_i] = res_i;
		}
	}
	
	return residuals;
}

void
MalamuteResidualOutput::outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols)
{
	
	std::string filename = getFilename();

	// Create file
	std::ofstream csv_file;
	csv_file.open(filename);
	
	// Column Names
	for (auto col_i : column_names){
		csv_file << col_i;
		if (col_i != column_names.back()){
			csv_file << ",";
		}
	}
	csv_file << " \n";
	
	// Body
	for (auto row_i : make_range(num_rows))
	{
		for (auto col_i : make_range(num_cols))
		{
			std::stringstream val;
			val << arr_2d[row_i][col_i];
			csv_file << val.str();
			csv_file << ",";
		}
		csv_file << " \n";
	}
	
	// Close file
	csv_file.close();
	
}

std::vector<Real>
MalamuteResidualOutput::get_x_y_z_id(const Point & p, const Real & id)
{
	std::vector<Real> x_y_z_id;
	x_y_z_id.push_back(p(0));  // x
	x_y_z_id.push_back(p(1));  // y
	x_y_z_id.push_back(p(2));  // z
	x_y_z_id.push_back(id);    // id
	
	return x_y_z_id;
	
}

std::string
MalamuteResidualOutput::getFilename()
{
	std::string filename = _file_base + "_" + _output_variable;
	std::string suffix;
	std::stringstream suffix_ss;
	
	int width = 4;
	
	if (_current_execute_flag == "TIMESTEP_BEGIN") {
		suffix_ss << std::setw(width) << std::setfill('0') << _t_step - 1;
		suffix = suffix_ss.str() + "_TIMESTEP_BEGIN";
	} else if (_current_execute_flag == "TIMESTEP_END") {
		suffix_ss << std::setw(width) << std::setfill('0') << _t_step - 1;
		suffix = suffix_ss.str() + "_TIMESTEP_END";
	} else if (_current_execute_flag == "NONLINEAR") {
		suffix_ss << std::setw(width) << std::setfill('0') << _t_step - 1;
		suffix_ss << "_";
		suffix_ss << std::setw(width) << std::setfill('0') << _nonlinear_iter;
		suffix = suffix_ss.str();
	}
	
	filename += "_" + suffix + ".csv";
	
	return filename;
		

}