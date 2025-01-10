#include "FVSystemOutputBase.h"

registerMooseObjectAliased("BubbleFlowApp", FVSystemOutputBase, "FVSystemOutputBase");

InputParameters
FVSystemOutputBase::validParams()
{
    InputParameters params = FVOutputBase::validParams();

    params.addParam<std::string>("which_time", "none", "Which times to display {none, all, converged, nonlinear, linear}");
	params.addParam<std::string>("time_units", "s", "Display-time units {s, ms, seconds, milliseconds}");

    return params;
}

FVSystemOutputBase::FVSystemOutputBase(const InputParameters & parameters)
  : FVOutputBase(parameters),
    dof_object(parameters),
    // Timing Inputs
    _which_time(getParam<std::string>("which_time")),
	_time_units(getParam<std::string>("time_units"))
{
    // Setup
    setupTimeVariables();
    setVariablesList();
}

void
FVSystemOutputBase::output()
{}

void
FVSystemOutputBase::setupTimeVariables()
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

void
FVSystemOutputBase::setVariablesList()
{
    if (_output_variable_name != "all")
        varaible_names.push_back(_output_variable_name);
    else
        for (auto v : _nl.getVariableNames())
            varaible_names.push_back(v);
}

void
FVSystemOutputBase::setDofMaps()
{
    for (std::string variable : varaible_names)
    {
        dof_maps[variable] = dof_object.getSingleVariableDOFMap(variable);
        dof_column_names[variable] = dof_object.getDOFMapColumnNames(dof_maps[variable]);
    }
}

FVDofTupleMap
FVSystemOutputBase::extractDofTuples(std::string variable_name, std::string component_name)
{
    // Get the DOF Map
    FVDofMap dof_map = dof_maps[variable_name];
    std::vector<std::string> column_names = dof_column_names[variable_name];

    // Get the possible DOF components {x, y, z} or {NONE (scalar)}
    std::pair<std::vector<std::string>, std::vector<std::string>> dof_var_comps = 
                                extractPossibleVectorComponents(column_names, variable_name);
    std::vector<std::string> dof_comps_full = dof_var_comps.first;
    std::vector<std::string> variable_comps_full = dof_var_comps.second;

    // Get the DOF components to be used
    std::vector<std::string> dof_comps;
    std::vector<std::string> var_comps;
    if (component_name == "all") {
        dof_comps = dof_comps_full;
        var_comps = variable_comps_full;
    }
    else {
        if (isInVector(dof_comps_full, "dof_" + component_name)) {
            dof_comps = {"dof_" + component_name}; 
            var_comps = {variable_name + "_" + component_name};
        }
        else {
            std::string error_message = "Error: dof component " + component_name + " not in { ";
            for (auto dof_c : dof_comps_full) error_message += dof_c + " "; error_message += "}";
            mooseError(error_message);
        }
    }

    // Combine to new map
    FVDofTupleMap dof_tuples;
    for (auto it=dof_map.begin(); it!=dof_map.end(); it++) {
        int elem_idx = it->first;
        std::map<std::string, Real> dofs_n = it->second;
        for (std::string dof_c : dof_comps) {
            int dof_i = dofs_n[dof_c];
            dof_tuples[dof_i] = std::make_tuple(elem_idx, dof_c);
        }
    }

    return dof_tuples;
}

std::pair<std::vector<std::string>, std::vector<std::string>>
FVSystemOutputBase::extractPossibleVectorComponents(std::vector<std::string> column_names, std::string variable_name)
{
    // Get the possible DOF components {x, y, z} or {NONE (scalar)}
    int number_of_geometric_columns = 4;
    std::vector<std::string> vec_comps = {"x", "y", "z"};
    
    int Ncomps = column_names.size() - number_of_geometric_columns;
    
    std::vector<std::string> variable_comps;
    std::vector<std::string> dof_comps;
    
    if (Ncomps == 1) {
        variable_comps.push_back(variable_name);
        dof_comps.push_back("dof");
    }
    else {
        for (int i=0; i<Ncomps; i++) {
            variable_comps.push_back(variable_name + "_" + vec_comps[i]);
            dof_comps.push_back("dof_" + vec_comps[i]);
        }
    }

    return std::make_pair(dof_comps, variable_comps);
}
