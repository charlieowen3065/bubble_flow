#include "LinearSystemBaseOutput.h"

registerMooseObjectAliased("BubbleFlowApp", LinearSystemBaseOutput, "LinearSystemBaseOutput");

InputParameters
LinearSystemBaseOutput::validParams()
{
    InputParameters params = FileOutput::validParams();

    params.addParam<NonlinearSystemName>("nl_sys", "nl0", "The nonlinear system that we should output information for.");
	params.addParam<std::string>("system_name", "nl0", "System to output");
	
	params.addParam<std::string>("variable", "all", "Name of output variable");

    params.addParam<bool>("output_DOF_map", true, "Boolean to output the DOF mapping");

    params.addParam<std::string>("which_time", "none", "Which times to display {none, all, converged, nonlinear, linear}");
	params.addParam<std::string>("time_units", "s", "Display-time units {s, ms, seconds, milliseconds}");

    return params;
}

LinearSystemBaseOutput::LinearSystemBaseOutput(const InputParameters & parameters)
  : FileOutput(parameters),
    // Nonlinear system variables
    _tid(getParam<THREAD_ID>("_tid")),
    _nl_sys_num(_problem_ptr->nlSysNum(getParam<NonlinearSystemName>("nl_sys"))),
	_nl(_problem_ptr->getNonlinearSystemBase(_nl_sys_num)),
    // Nonlinear-Variable Information
	_output_variable_name(getParam<std::string>("variable")),
    // System variables
    _sys(_problem_ptr->es().get_system(getParam<std::string>("system_name"))),
	_dof_map(_sys.get_dof_map()),
	_sys_number(_dof_map.sys_number()),
	_system_name(getParam<std::string>("system_name")),
    // Meshing
    _mesh(_problem_ptr->mesh()),
	_num_nodes(_mesh.nNodes()),
    // Timing Inputs
    _which_time(getParam<std::string>("which_time")),
	_time_units(getParam<std::string>("time_units"))
{
    if (_output_variable_name != "all")
    {
        if (!(_problem_ptr->hasVariable(_output_variable_name)))
            mooseError("ERROR: Problem does not contain the input variable | " + _output_variable_name);
    }
    // Compiles booleans to dispaly times
    setupTimeVariables();
	
}

void
LinearSystemBaseOutput::output() {}

void
LinearSystemBaseOutput::setupTimeVariables()
{
    if (_which_time == "none"){
		_show_converged_time = false;
		_show_nl_time = false;
		_show_l_time = false;
	} else if (_which_time == "all"){
		_show_converged_time = true;
		_show_nl_time = true;
		_show_l_time = true;
	} else if (_which_time == "converged"){
		_show_converged_time = true;
		_show_nl_time = false;
		_show_l_time = false;
	} else if (_which_time == "nonlinear"){
		_show_converged_time = true;
		_show_nl_time = true;
		_show_l_time = false;
	} else if (_which_time == "linear"){
		_show_converged_time = true;
		_show_nl_time = true;
		_show_l_time = true;
	} 
}

std::string
LinearSystemBaseOutput::getFilename(std::string BASE)
{
    std::string filename = BASE + "_" + _output_variable_name;

    int width = 4;
    std::stringstream suffix_ss;

    suffix_ss << "t_" << std::setw(width) << std::setfill('0') << _t_step;  // time-step
    if (_current_execute_flag == "TIMESTEP_BEGIN") {
        suffix_ss << "_TIMESTEP_BEGIN";
    } else if (_current_execute_flag == "TIMESTEP_END") {
        suffix_ss << "_TIMESTEP_END";
    } else {
        suffix_ss << "_";
        suffix_ss << "nl_" << std::setw(width) << std::setfill('0') << _nonlinear_iter;   // Nonlinear iteration
		suffix_ss << "_";
		suffix_ss << "l_" << std::setw(width) << std::setfill('0') << _linear_iter;       // Linear iteration

        if (_current_execute_flag == "NONLINEAR") {
            suffix_ss << "_NONLINEAR";
        } else if (_current_execute_flag == "LINEAR") {
            suffix_ss << "LINEAR";
        }
    }
    std::string suffix = suffix_ss.str();
    filename += "_" + suffix;
    return filename;
}

void
LinearSystemBaseOutput::outputToCSV(Real** arr_2d, std::vector<std::string> column_names, 
                                    int num_rows, int num_cols, std::string base_name)
{
	
	std::string filename = getFilename(base_name) + ".csv";

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
LinearSystemBaseOutput::getNodeCoordinates(const Point & p, const Real & id)
{
	std::vector<Real> coords;
	coords.push_back(p(0));  // x
	coords.push_back(p(1));  // y
	coords.push_back(p(2));  // z
	coords.push_back(id);    // id
	return coords;
}

dof_id_type**
LinearSystemBaseOutput::getDOFIndices_SingleVariable(std::vector<dof_id_type> & di, const unsigned int sys_number,
                                           			 const unsigned int var_number, unsigned int number_componets)
{
	// Declare return array
	dof_id_type** component_array = new dof_id_type*[number_componets];
	// Get dof indices
	for (auto comp_number : make_range(number_componets))
	{
		component_array[comp_number] = new dof_id_type[_num_nodes];
		for (auto node_idx : make_range(_num_nodes))
		{
			const Node * node_i = _mesh.nodePtr(node_idx);
			dof_id_type dof_i = node_i->dof_number(sys_number, var_number, comp_number);
			di.push_back(dof_i);
			component_array[comp_number][node_idx] = dof_i;		
		}
	}
	return component_array;
}