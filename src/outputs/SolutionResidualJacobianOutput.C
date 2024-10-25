#include "SolutionResidualJacobianOutput.h"

registerMooseObjectAliased("MalamuteApp", SolutionResidualJacobianOutput, "SolutionResidualJacobianOutput");

InputParameters
SolutionResidualJacobianOutput::validParams()
{
	InputParameters params = FileOutput::validParams();
	
	params.addParam<NonlinearSystemName>("nl_sys", "nl0", "The nonlinear system that we should output information for.");
	params.addParam<std::string>("system_name", "nl0", "System to output");
	
	params.addParam<std::string>("variable", "Name of output variable");
	
	// Solution, Residual, Jacobian Inputs
	params.addParam<bool>("output_solution", true, "Boolean to output the solution vector");
	params.addParam<bool>("output_delta_solution", true, "Boolean to output the change in solution vector");
	params.addParam<bool>("output_resdiual", true, "Boolean to output the resdiual vector");
	params.addParam<bool>("output_jacobian", true, "Boolean to output the jacobian vector");
	params.addParam<bool>("output_DOF_map", true, "Boolean to output the DOF mapping");
	
	// Timing Inputs
	params.addParam<std::string>("which_time", "none", "Which times to display {none, all, converged, nonlinear, linear}");
	params.addParam<std::string>("time_units", "s", "Display-time units {s, ms, seconds, milliseconds}");
	
	// Solution & Residual Specific Booleans
	params.addParam<bool>("seperate_solution_and_residual", false, "creates a file for both the solution and residual instead of just 1 file");
	
	// Jacobian-Specific Inputs
	params.addParam<bool>("display_jacobian", false, "Boolean to control if the Jacobian is printed to the screen");
	params.addParam<bool>("save_jacobian_txt", true, "Boolean to control if the Jacobian is saved to a txt file");
	params.addParam<bool>("save_jacobian_csv", true, "Boolean to control if the Jacobian is saved to a csv file");
	
	// Debug
	params.addParam<bool>("debug_SolutionResidualJacobianOutput", false, "Turns debug on and off");
	params.addParam<bool>("show_L2_norm", false, "DEBUG: Show the L2-Norm of the residual vector");

	
	return params;
}

SolutionResidualJacobianOutput::SolutionResidualJacobianOutput(const InputParameters & parameters)
  : FileOutput(parameters),
    
	// Nonlinear System Variables
    _tid(getParam<THREAD_ID>("_tid")),
	_nl_sys_num(_problem_ptr->nlSysNum(getParam<NonlinearSystemName>("nl_sys"))),
	_nl(_problem_ptr->getNonlinearSystemBase(_nl_sys_num)),
	
	// System Variables
	_sys(_problem_ptr->es().get_system(getParam<std::string>("system_name"))),
	_dof_map(_sys.get_dof_map()),
	_sys_number(_dof_map.sys_number()),
	_system_name(getParam<std::string>("system_name")),
	
	// Meshing
	_mesh(_problem_ptr->mesh()),
	_num_nodes(_mesh.nNodes()),
	
	// Nonlinear-Variable Information
	_output_variable_name(getParam<std::string>("variable")),
	_output_variable(_nl.getVariable(_tid, _output_variable_name)),
	
	// Solution, Residual, and Jacobian References
	_solution(_nl.solution()),
	_residual(_nl.RHS()),
	residual_full(static_cast<ImplicitSystem &>(_nl.system()).get_vector(0)),
	
	_prev_it_solution(0),
	_prev_it_residual(0),
	
	_linear_itt_counter(0),
	_copy_var(true),
	
	// ---------- User-Inputs ---------- //
	// Solution, Residual, Jacobian Inputs
	_output_solution(getParam<bool>("output_solution")),
	_output_delta_solution(getParam<bool>("output_delta_solution")),
	_output_residual(getParam<bool>("output_resdiual")),
	_output_jacobian(getParam<bool>("output_jacobian")),
	_output_DOF_map(getParam<bool>("output_DOF_map")),
	// Timing Inputs
	_which_time(getParam<std::string>("which_time")),
	_time_units(getParam<std::string>("time_units")),
	// Solution & Residual Specific Inputs
	_solution_residual_seperate_files(getParam<bool>("seperate_solution_and_residual")),
	// Jacobian-Specific Inputs
	_display_jacobian(getParam<bool>("display_jacobian")),
	_save_jacobian_txt(getParam<bool>("save_jacobian_txt")),
	_save_jacobian_csv(getParam<bool>("save_jacobian_csv")),
	
	// Debug
	_SRJ_debug(getParam<bool>("debug_SolutionResidualJacobianOutput")),
	_show_L2_norm(getParam<bool>("show_L2_norm"))
	
{	
	// Checks that the problem contains the desired variable
	if (!(_problem_ptr->hasVariable(_output_variable_name)))
	{
		mooseError("ERROR: Problem does not contain the input variable | " + _output_variable_name);
	}
	
	// Compiles booleans to dispaly times
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

// ================================================================ MAIN OUTPUT =============================================================== //

void
SolutionResidualJacobianOutput::output()
{	
	_output_variable_number = _output_variable.number();
	_number_componets = _mesh.nodePtr(0)->n_comp(_sys_number, _output_variable_number);
	if (_number_componets == 1) {_is_vector = false;} else {_is_vector = true;}
	
	displayTime();
	
	if (_output_DOF_map)
		outputDOFMap();
	
	// --------- Solution & Residual --------- //
	outputSolutionResidual();
	
	if (_output_delta_solution)
		outputDeltaSolution();
	
	// -------------- Jacobian --------------- //
	if (_output_jacobian)
	{
		JacobianOutput();
	}
	
	if (_SRJ_debug)
	{
		debugOutput();
	}
	
}

// =============================================================== MISC. METHODS =============================================================== //


std::string
SolutionResidualJacobianOutput::getFilename(std::string BASE)
{
	std::string filename = BASE + "_" + _output_variable_name;
	
	int width = 4;
	std::stringstream suffix_ss;
		
	suffix_ss << "t_" << std::setw(width) << std::setfill('0') << _t_step;        // Time-step
	if (_current_execute_flag == "TIMESTEP_BEGIN") {
		suffix_ss << "_TIMESTEP_BEGIN";
	} else if (_current_execute_flag == "TIMESTEP_END") {
		suffix_ss << "_TIMESTEP_END";
	} else {
		suffix_ss << "_";
		suffix_ss << "nl_" << std::setw(width) << std::setfill('0') << _nonlinear_iter;   // Nonlinear iteration
		suffix_ss << "_";
		suffix_ss << "l_" << std::setw(width) << std::setfill('0') << _linear_iter;       // Linear iteration
		
		if (_current_execute_flag == "NONLINEAR"){
			suffix_ss << "_NONLINEAR";
		} else if (_current_execute_flag == "LINEAR"){
			suffix_ss << "_LINEAR";
		}
		
	}
	std::string suffix = suffix_ss.str();
	
	
	
	filename += "_" + suffix;
	
	return filename;
}

dof_id_type**
SolutionResidualJacobianOutput::getDofIndices(std::vector<dof_id_type> & di, const unsigned int sys_number, const unsigned int var_number)
{
	
	// Declare return array
	dof_id_type** component_array = new dof_id_type*[_number_componets];
	
	
	// Get dof indices
	for (auto comp_number : make_range(_number_componets))
	{
		component_array[comp_number] = new dof_id_type[_num_nodes];
		
		for (auto node_idx : make_range(_num_nodes))
		{
			const Node * node_i = _mesh.nodePtr(node_idx);
			
			dof_id_type dof_i = node_i->dof_number(_sys_number, _output_variable_number, comp_number);
			
			di.push_back(dof_i);
			component_array[comp_number][node_idx] = dof_i;
			
		}
		
	}
	
	return component_array;
}

void
SolutionResidualJacobianOutput::getDOFMapping(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	
	std::vector<dof_id_type> dof_indices_full;
	dof_id_type** component_array = getDofIndices(dof_indices_full, _sys_number, _output_variable_number);
	
	// Place DOF Mapping into vector
	if (_is_vector == false) {
		if (_output_solution){
			variable_names.push_back("Solution");
		}
		if (_output_residual){
			variable_names.push_back("Residual");
		}
		std::vector<dof_id_type> dof_map_i;
		for (auto i : dof_indices_full) {
			dof_map_i.push_back(i);
		}
		
		dof_mapping.push_back(dof_map_i);	
	} else {
		
		std::vector<std::string> vec_comps = {"_x", "_y", "_z"};
		
		// Variable Names -- Solution
		if (_output_solution){
			for (auto c_i : make_range(_number_componets))
			{
				variable_names.push_back("Solution" + vec_comps[c_i]);		
			}
		}
		// Variable Names -- Residual
		if (_output_residual){
			for (auto c_i : make_range(_number_componets))
			{
				variable_names.push_back("Residual" + vec_comps[c_i]);	
			}
		}
		// DOF Mapping
		for (auto c_i : make_range(_number_componets))
		{
			std::vector<dof_id_type> dof_map_i;
			for (auto i : make_range(_num_nodes)) 
			{
				dof_map_i.push_back(component_array[c_i][i]);
			}
		
			dof_mapping.push_back(dof_map_i);		
		}	
	}

}

void
SolutionResidualJacobianOutput::outputDOFMap()
{
	std::vector<std::string> variable_names;
	std::vector<std::vector<dof_id_type>> dof_mapping;
	
	getDOFMapping(variable_names, dof_mapping);
	
	std::string base_filename = "DOF_Mapping";
	outputToCSV(dof_mapping, variable_names, _num_nodes, _number_componets, base_filename);
}

void
SolutionResidualJacobianOutput::displayTime()
{
	// Initial Checks
	
	if (_which_time == "none"){return;}
	if ((_which_time != "linear") and (_which_time != "all"))
	{
		if (_current_execute_flag == "LINEAR") {return;}
		if (_which_time != "nonlinear")
		{
			if (_current_execute_flag == "NONLINEAR") {return;}
		}
	}
	
	// Time-Units
	std::string units;
	if ((_time_units == "s") or (_time_units == "seconds"))
	{
		units = "s";
	} else if ((_time_units == "ms") or (_time_units == "milliseconds"))
	{
		units = "ms";
	}
	
	// Exce. Flag Name
	std::string flag_name;
	if (_current_execute_flag == "TIMESTEP_BEGIN")
	{
		flag_name = "T-Step Begin";
	} else if (_current_execute_flag == "TIMESTEP_END")
	{
		flag_name = "T-Step End";
	} else if (_current_execute_flag == "NONLINEAR")
	{
		flag_name = "Nonlinear";
	} else if (_current_execute_flag == "LINEAR")
	{
		flag_name = "Linear";
	}
		
	// Display time and flag
	int NOW_TIME_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	int time_ms = NOW_TIME_ms - START_TIME_ms;
	
	if (units == "s")
	{
		double time_s = (double)time_ms / 1000.0;
		std::cout << "*** " << flag_name << " Time: " << time_s << " [s] *** " << std::endl;
	} else if (units == "ms")
	{
		std::cout << "*** " << flag_name << " Time: " << time_ms << " [ms] *** " << std::endl;		
	}
}

Real
SolutionResidualJacobianOutput::computeL2Norm(std::vector<Real> residual_vector)
{
	Real squared_sum = 0;
	for (Real r_i : residual_vector)
		squared_sum += (r_i * r_i);
	Real l2_norm = 0.5 * std::pow(squared_sum, 0.5);
	return l2_norm;
}


// ================================================ OUTPUT METHODS FOR THE SOLUTION & RESIDUALS ================================================ //
void
SolutionResidualJacobianOutput::outputSolutionResidual()
{
	// --- Solution & Residuals --- //
	// 4 Cases:
	//		1) Soln & Res (together)
	//		2) Soln & Res (seperate)
	//		3) Soln only
	//		4) Res only
	
	computeLinearIterationOutputs();
	
	std::vector<std::string> variable_names;
	std::vector<std::vector<dof_id_type>> dof_mapping;
	Real** SolutionResidual;
	std::vector<std::string> column_names;
	std::string base_filename;
	
	if ((_output_solution or _output_residual) and (!_solution_residual_seperate_files)){
		getDOFMapping(variable_names, dof_mapping);
		
		SolutionResidual = createSolutionResidualOutputArray(variable_names, dof_mapping);
		base_filename = "SolutionResidual";
		
		column_names.push_back("id");
		column_names.push_back("x");
		column_names.push_back("y");
		column_names.push_back("z");
		for (auto out_name : variable_names) {
			column_names.push_back(out_name);
		}
		
		int addon1;
		if (_output_solution and _output_residual){
			addon1 = 2 * _number_componets;
		} else if (_output_solution){
			addon1 = _number_componets;
		} else if (_output_residual){
			addon1 = _number_componets;
		}
		
		outputToCSV(SolutionResidual, column_names, _num_nodes, 4 + addon1, base_filename);
		
	} else if ((_output_solution and _output_residual) and (_solution_residual_seperate_files)){
		// Solution only
		_output_solution = true;
		_output_residual = false;
		
		getDOFMapping(variable_names, dof_mapping);
		
		SolutionResidual = createSolutionResidualOutputArray(variable_names, dof_mapping);
		base_filename = "Solution";
		
		column_names.push_back("id");
		column_names.push_back("x");
		column_names.push_back("y");
		column_names.push_back("z");
		for (auto out_name : variable_names) {
			column_names.push_back(out_name);
		}
		
		outputToCSV(SolutionResidual, column_names, _num_nodes, 4 + _number_componets, base_filename);
		
		// Residual only
		_output_solution = false;
		_output_residual = true;
		variable_names.clear();
		column_names.clear();
		dof_mapping.clear();
		
		getDOFMapping(variable_names, dof_mapping);
		
		SolutionResidual = createSolutionResidualOutputArray(variable_names, dof_mapping);
		base_filename = "Residual";
		
		column_names.push_back("id");
		column_names.push_back("x");
		column_names.push_back("y");
		column_names.push_back("z");
		for (auto out_name : variable_names) {
			column_names.push_back(out_name);
		}
		
		outputToCSV(SolutionResidual, column_names, _num_nodes, 4 + _number_componets, base_filename);
		
		_output_solution = true;
		_output_residual = true;
		
	}
	
	else if (_output_solution){
		// Solution only
		_output_solution = true;
		_output_residual = false;
		
		getDOFMapping(variable_names, dof_mapping);
		
		SolutionResidual = createSolutionResidualOutputArray(variable_names, dof_mapping);
		base_filename = "Solution";
		
		column_names.push_back("id");
		column_names.push_back("x");
		column_names.push_back("y");
		column_names.push_back("z");
		for (auto out_name : variable_names) {
			column_names.push_back(out_name);
		}
		
		outputToCSV(SolutionResidual, column_names, _num_nodes, 4 + _number_componets, base_filename);
		
		_output_solution = true;
		
	}
}

void
SolutionResidualJacobianOutput::outputDeltaSolution()
{
	//std::cout << "Here 1" << std::endl;
	
	std::vector<std::string> variable_names;
	std::vector<std::vector<dof_id_type>> dof_mapping;
	Real** DeltaSolution;
	std::vector<std::string> column_names;
	std::string base_filename = "DeltaSolution";
	
	getDOFMapping(variable_names, dof_mapping);
	DeltaSolution = createDeltaSolutionOutputArray(variable_names, dof_mapping);
	
	column_names.push_back("id");
	column_names.push_back("x");
	column_names.push_back("y");
	column_names.push_back("z");
	
	for (auto out_name : variable_names) {
		column_names.push_back(out_name);
	}
	
	int addon = _number_componets;
	
	outputToCSV(DeltaSolution, column_names, _num_nodes, 4 + addon, base_filename);
}

Real**
SolutionResidualJacobianOutput::createSolutionResidualOutputArray(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	Real** OutputArray = new Real*[_num_nodes];
	
	const NumericVector<Real> * curr_solution = _nl.currentSolution();
	
	int addon1;
	int addon2;
	
	if (_output_solution and _output_residual){
		addon1 = 2 * _number_componets;
		addon2 = _number_componets;
	} else if (_output_solution){
		addon1 = _number_componets;
		addon2 = 0;
	} else if (_output_residual){
		addon1 = _number_componets;
		addon2 = 0;
	}

	
	for (dof_id_type node_idx : make_range(_num_nodes))
	{
		const Node * node_i = _mesh.nodePtr(node_idx);
		std::vector<Real> x_y_z_id = get_x_y_z_id(*node_i, node_idx);
		
		OutputArray[node_idx] = new Real[4 + addon1];
		OutputArray[node_idx][0] = x_y_z_id[3];  // id
		OutputArray[node_idx][1] = x_y_z_id[0];  // x
		OutputArray[node_idx][2] = x_y_z_id[1];  // y
		OutputArray[node_idx][3] = x_y_z_id[2];  // z
		
		// Solution
		if (_output_solution){
			for (auto c_i : make_range(_number_componets))
			{
				dof_id_type dof_i = dof_mapping[c_i][node_idx];
				//Real soln_i = _solution(dof_i);
				Real soln_i = curr_solution->el(dof_i);
				
				OutputArray[node_idx][4 + c_i] = soln_i;
			}
		}
		// Residual
		if (_output_residual){
			for (auto c_i : make_range(_number_componets))
			{
				dof_id_type dof_i = dof_mapping[c_i][node_idx];
				Real res_i = _residual(dof_i);
				
				OutputArray[node_idx][4 + addon2 + c_i] = res_i;
			}
		}
	}
	
	return OutputArray;
}

Real**
SolutionResidualJacobianOutput::createDeltaSolutionOutputArray(std::vector<std::string> & variable_names, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	//std::cout << "Here 2" << std::endl;
	
	Real** OutputArray = new Real*[_num_nodes];
	
	int addon = _number_componets;
	NumericVector<Real> * prev_solution = _nl.solutionPreviousNewton();
	
	for (dof_id_type node_idx : make_range(_num_nodes))
	{
		const Node * node_i = _mesh.nodePtr(node_idx);
		std::vector<Real> x_y_z_id = get_x_y_z_id(*node_i, node_idx);
		
		OutputArray[node_idx] = new Real[4 + addon];
		OutputArray[node_idx][0] = x_y_z_id[3];  // id
		OutputArray[node_idx][1] = x_y_z_id[0];  // x
		OutputArray[node_idx][2] = x_y_z_id[1];  // y
		OutputArray[node_idx][3] = x_y_z_id[2];  // z
		
		// Delta-Solution
		if (_output_solution){
			for (auto c_i : make_range(_number_componets))
			{
				std::cout << "Here 21" << std::endl;
				dof_id_type dof_i = dof_mapping[c_i][node_idx];
				
				std::cout << "Here 22 (" << dof_i << ")" << std::endl;
				Real soln_i = _solution(dof_i);
				std::cout << "prev_solution->size(): " << prev_solution->size() << std::endl;
				Real prev_soln_i = prev_solution->el(dof_i);
				std::cout << "Here 23" << std::endl;
				Real delta_S_i = soln_i - prev_soln_i;
				std::cout << "Here 24" << std::endl;

				OutputArray[node_idx][4 + c_i] = delta_S_i;
				std::cout << "Here 25" << std::endl;
			}
		}
	}
	
	return OutputArray;
}

std::vector<Real>
SolutionResidualJacobianOutput::getResidualVector_SingleComponent(int component_number, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	std::vector<Real> residual_vector;
	_problem_ptr->computeResidual(_solution, residual_full, _nl_sys_num);
	
	
	for (dof_id_type node_idx : make_range(_num_nodes))
	{
		dof_id_type dof_i = dof_mapping[component_number][node_idx];
		Real res_i = residual_full(dof_i);
		residual_vector.push_back(res_i);
	}
	
	return residual_vector;
}

std::vector<Real>
SolutionResidualJacobianOutput::getSolutionVector_SingleComponent(int component_number, std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	std::vector<Real> solution_vector;	
	
	for (dof_id_type node_idx : make_range(_num_nodes))
	{
		dof_id_type dof_i = dof_mapping[component_number][node_idx];
		Real soln_i = _solution(dof_i);
		solution_vector.push_back(soln_i);
	}
	
	return solution_vector;
}

void
SolutionResidualJacobianOutput::computeLinearIterationOutputs()
{	
	if (_current_execute_flag == "LINEAR")
		mooseError("ERROR: Cannot do linear iterations at the moment.");
	
	if (_current_execute_flag == "LINEAR")
	{
		_linear_itt_counter += 1;
		if (_copy_var == true)
		{
			getPreviousVectors();
			_copy_var = false;
		}
	} 
	else {
		// Booleans on if to copy the vectors next Linear-step
		if (_current_execute_flag == "NONLINEAR")
			_copy_var = true;
		
		if (_linear_itt_counter > 0)
		{
			allLinearIterations(_linear_itt_counter);
			_linear_itt_counter = 0;
		}
	}
	
}

void
SolutionResidualJacobianOutput::allLinearIterations(int number_linear_iterations)
{
	std::cout << "\n\nallLinearIterations" << std::endl;
	std::cout << "Num. Linear: " << number_linear_iterations << std::endl;
	
	Real l2_norm = computeL2Norm(_prev_it_residual);
	std::cout << "l2_norm: " << l2_norm * 2 << "\n" << std::endl;
	
	//NumericVector<Real> & residual_use = static_cast<NumericVector<Real>>(_prev_it_residual);
	//NumericVector<Real> & solution_use = static_cast<NumericVector<Real>>(_prev_it_solution);
	
	
		
}

void
SolutionResidualJacobianOutput::getPreviousVectors()
{	
	_prev_it_residual.clear();
	_prev_it_solution.clear();
	
	for (int i=0; i < _residual.size(); i++)
		_prev_it_residual.push_back(_residual.el(i));
	for (int i=0; i < _solution.size(); i++)
		_prev_it_solution.push_back(_solution.el(i));
}


void
SolutionResidualJacobianOutput::outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name)
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

void
SolutionResidualJacobianOutput::outputToCSV(std::vector<std::vector<dof_id_type>> arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name)
{
	std::cout << "START" << std::endl;
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
			//std::cout << row_i << ", " << col_i << " | val: " << arr_2d[row_i][col_i] << std::endl;
			std::stringstream val;
			val << arr_2d[col_i][row_i];
			csv_file << val.str();
			csv_file << ",";
		}
		csv_file << " \n";
	}
	
	// Close file
	csv_file.close();
	std::cout << "END" << std::endl;
	
}

std::vector<Real>
SolutionResidualJacobianOutput::get_x_y_z_id(const Point & p, const Real & id)
{
	std::vector<Real> x_y_z_id;
	x_y_z_id.push_back(p(0));  // x
	x_y_z_id.push_back(p(1));  // y
	x_y_z_id.push_back(p(2));  // z
	x_y_z_id.push_back(id);    // id
	
	return x_y_z_id;
	
}

// ====================================================== OUTPUT METHODS FOR THE JACOBIAN ====================================================== //
void
SolutionResidualJacobianOutput::JacobianOutput()
{
	// Assemble matrix
	// Petsc stuff
	auto & petsc_options = _problem_ptr->getPetscOptions();
	auto & pars = _problem_ptr->solverParams();
	Moose::PetscSupport::petscSetOptions(petsc_options, pars);
	ExecFlagType flag = _current_execute_flag;
	_problem_ptr->execute(flag);
	
	auto & jacobian = static_cast<ImplicitSystem &>(_nl.system()).get_system_matrix();
	
	//_problem_ptr->computeJacobian(*_nl.currentSolution(), _jacobian, _nl_sys_num);
	_problem_ptr->computeJacobian(*_nl.currentSolution(), jacobian, _nl_sys_num);
	
	if (_save_jacobian_txt)
	{
		std::cout << "SAVING JACOBIAN" << std::endl;
		saveJacobianToTextFile(jacobian);
	}
	
	if (_display_jacobian)
	{
		std::cout << "JACOBIAN INFORMATION: " << std::endl;
		jacobian.print(std::cout, true);
	}
	
}

void
SolutionResidualJacobianOutput::saveJacobianToTextFile(SparseMatrix<Real> & jacobian)
{
	std::string base_name = "Jacobian";
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
SolutionResidualJacobianOutput::saveJacobianToCSVFile(SparseMatrix<Real> & jacobian)
{
	std::string base_name = "Jacobian";
	std::string filename = getFilename(base_name) + ".csv";
	
	// Create file
	std::ofstream csv_file;
	csv_file.open(filename);
	
}

// =============================================================== DEBUG METHODS =============================================================== //

void
SolutionResidualJacobianOutput::debugOutput()
{
	std::vector<std::string> variable_names;
	std::vector<std::vector<dof_id_type>> dof_mapping;
	getDOFMapping(variable_names, dof_mapping);
	
	if (_show_L2_norm)
	{
		for (auto c_i : make_range(_number_componets))
		{
			std::string var_name = variable_names[c_i];
			std::vector<Real> res_vec = getResidualVector_SingleComponent(c_i, dof_mapping);
			Real l2_norm = computeL2Norm(res_vec);
			//std::vector<Real> soln_vec = getSolutionVector_SingleComponent(c_i, dof_mapping);
			//Real l2_norm = computeL2Norm(soln_vec);
			std::cout << var_name << " L2-Norm: " << l2_norm << std::endl;
		}
	}
}









