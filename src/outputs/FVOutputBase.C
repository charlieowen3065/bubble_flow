#include "FVOutputBase.h"

registerMooseObjectAliased("BubbleFlowApp", FVOutputBase, "FVOutputBase");

InputParameters
FVOutputBase::validParams()
{
    InputParameters params = FileOutput::validParams();

    params.addParam<NonlinearSystemName>("nl_sys", "nl0", "The nonlinear system that we should output information for.");
	params.addParam<std::string>("system_name", "nl0", "System to output");
	params.addParam<std::string>("variable", "all", "Name of output variable");
    params.addParam<bool>("output_DOF_map", false, "Boolean to output the DOF mapping");

    return params;
}

FVOutputBase::FVOutputBase(const InputParameters & parameters)
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
    _num_elems(_mesh.nElem())
{
    if (_output_variable_name != "all")
        if (!(_problem_ptr->hasVariable(_output_variable_name)))
            mooseError("ERROR: Problem does not contain the input variable | " + _output_variable_name);
}

void
FVOutputBase::output() 
{}

std::string
FVOutputBase::getFilename(std::string BASE)
{
    // std::string filename = BASE + "_" + _output_variable_name;
    std::string filename = BASE;

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
FVOutputBase::outputToCSV(Real** arr_2d, std::vector<std::string> column_names, 
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

void
FVOutputBase::vectorMapToCSV(FVVectorTupleMap vector_map, std::string base_filename)
{
    // Filename
    std::string filename = getFilename(base_filename);

    // Create file
	std::ofstream csv_file;
	csv_file.open(filename);

    // Column Names
    std::vector<std::string> column_names = {"Element", "DOF #", 
                                            "x", "y", "z", "id",
                                            "Variable Name", "Value"};
    for (auto col : column_names) {
        csv_file << col;
		if (col != column_names.back()){
			csv_file << ",";
		}
    } csv_file << " \n";

    // Extract DOF data
    for (auto it=vector_map.begin(); it!=vector_map.end(); it++){
        // Extract map data
        int dof_i     = it->first;
        auto dof_data = it->second;

        // Get DOF data
        int elem_idx         = std::get<0>(dof_data);
        std::string vec_comp = std::get<1>(dof_data);
        Real vector_value    = std::get<2>(dof_data);

        // Get {x, y, z, id} coordinates
        std::vector<Real> coords = getElementCoordinates(*_mesh.elemPtr(elem_idx), elem_idx);

        // Save data to CSV
        csv_file << elem_idx     << ",";  // Elem
        csv_file << dof_i        << ",";  // DOF #
        csv_file << coords[0]    << ",";  // x
        csv_file << coords[1]    << ",";  // y
        csv_file << coords[2]    << ",";  // z
        csv_file << coords[3]    << ",";  // id
        csv_file << vec_comp     << ",";  // variable name
        csv_file << vector_value;  // value
        csv_file << " \n";
    }

    // Close file
	csv_file.close();	
}

std::vector<Real>
FVOutputBase::getNodeCoordinates(const Point & p, const Real & id)
{
	std::vector<Real> coords;
	coords.push_back(p(0));  // x
	coords.push_back(p(1));  // y
	coords.push_back(p(2));  // z
	coords.push_back(id);    // id
	return coords;
}

std::vector<Real>
FVOutputBase::getElementCoordinates(const Elem & elem, const Real & id)
{
    Point p = elem.true_centroid();

	std::vector<Real> coords;
	coords.push_back(p(0));  // x
	coords.push_back(p(1));  // y
	coords.push_back(p(2));  // z
	coords.push_back(id);    // id
	return coords;
}

bool
FVOutputBase::isInVector(std::vector<std::string> vec, std::string target)
{
    return (find(vec.begin(), vec.end(), target) != vec.end());
}