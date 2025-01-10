#include "FVDOFMappingOutput.h"

registerMooseObjectAliased("BubbleFlowApp", FVDOFMappingOutput, "FVDOFMappingOutput");

InputParameters
FVDOFMappingOutput::validParams()
{
     InputParameters params = FVOutputBase::validParams();
     return params;
}

FVDOFMappingOutput::FVDOFMappingOutput(const InputParameters & parameters)
  : FVOutputBase(parameters)
{

}

void
FVDOFMappingOutput::output() 
{
    // Collects the variable names
    std::vector<VariableName> var_names;
    if ((_output_variable_name != "all") && (_output_variable_name != "none")) {
        var_names.push_back(_output_variable_name);
    } else if (_output_variable_name == "all") {
        for (auto v : _nl.getVariableNames())
            var_names.push_back(v);
    }

    // Outputs the DOF map for each given variable
    for (auto var_name : var_names){
        outputSingleDOFMap(var_name);
    }
}

dof_id_type**
FVDOFMappingOutput::getDOFIndices_SingleVariable(std::vector<dof_id_type> & di, const unsigned int sys_number,
                                           	     const unsigned int var_number, unsigned int number_componets)
{
    /* 
    Collects the dof indicies for an input variable, based on the variable number.

    Outputs a 2D array, where the first dimension is the components (this is 1D for
    a scalar variable) and the second dimension, the number of elements in the system.

    (NOTE: Elements instead of nodes as Finite-Volume variables are defined by elements)
    */
	
    // Declare return array
	dof_id_type** component_array = new dof_id_type*[number_componets];
	// Get dof indices
	for (auto comp_number : make_range(number_componets))
	{
		component_array[comp_number] = new dof_id_type[_num_elems];
		for (auto elem_idx : make_range(_num_elems))
		{
			const Elem * elem_i = _mesh.elemPtr(elem_idx);
			dof_id_type dof_i = elem_i->dof_number(sys_number, var_number, comp_number);
			di.push_back(dof_i);
			component_array[comp_number][elem_idx] = dof_i;		
		}
	}
	return component_array;
}

void
FVDOFMappingOutput::getSingleVariableDOFMapping(const unsigned int variable_number,
                                              std::vector<std::string> & dof_names, 
                                              std::vector<std::vector<dof_id_type>> & dof_mapping)
{
    // Determines if the variable is a vector or scalar
	bool is_vector;
    // unsigned int number_componets = _mesh.elemPtr(0)->n_comp(_sys_number, variable_number);
    unsigned int number_componets = 1; // For FV variables, all are scalars
    if (number_componets == 1) {is_vector = false;} else {is_vector = true;}

    // Set up & compute arrays of dofs
	std::vector<dof_id_type> dof_indices_full;
	dof_id_type** component_array = getDOFIndices_SingleVariable(dof_indices_full, _sys_number,
                                                                 variable_number, number_componets);

    if (!is_vector) // Not a vector variable -> only 1 set of dofs
    {
        dof_names.push_back("dof");
        std::vector<dof_id_type> dof_map_i;
        for (auto i : dof_indices_full)
            dof_map_i.push_back(i);
        dof_mapping.push_back(dof_map_i);
    } 
    else  // is a vector variable -> 1-3 sets of dofs ({x, y, z}-components has seperate dofs)
    {
        std::vector<std::string> vec_comps = {"x", "y", "z"};
        for (auto c_i : make_range(number_componets))
		{
            dof_names.push_back("dof_" + vec_comps[c_i]);
			std::vector<dof_id_type> dof_map_i;
			for (auto i : make_range(_num_elems)) 
				dof_map_i.push_back(component_array[c_i][i]);
			dof_mapping.push_back(dof_map_i);		
		}
    }
}

Real**
FVDOFMappingOutput::createDOFMapOutputArray(std::vector<std::string> & dof_names,
                                          std::vector<std::vector<dof_id_type>> & dof_mapping)
{
    Real** OutputArray = new Real*[_num_elems];
    int _number_componets = dof_names.size();

    for (dof_id_type elem_idx : make_range(_num_elems))
    {
        const Elem * elem_i = _mesh.elemPtr(elem_idx);
        std::vector<Real> coords = getElementCoordinates(*elem_i, elem_idx);

        OutputArray[elem_idx] = new Real[4 + _number_componets];
        OutputArray[elem_idx][0] = coords[0];
        OutputArray[elem_idx][1] = coords[1];
        OutputArray[elem_idx][2] = coords[2];
        OutputArray[elem_idx][3] = coords[3];
        
        for (auto ci : make_range(_number_componets))
        {
            dof_id_type dof_i = dof_mapping[ci][elem_idx];
            OutputArray[elem_idx][4 + ci] = dof_i;
        }
    }
    return OutputArray;
}

void
FVDOFMappingOutput::outputSingleDOFMap(std::string variable_name)
{
    int variable_number = _nl.getVariable(_tid, variable_name).number();
    std::vector<std::string> dof_names;
    std::vector<std::vector<dof_id_type>> dof_mapping;
    getSingleVariableDOFMapping(variable_number, dof_names, dof_mapping);

    Real** dofMapArray;
    dofMapArray = createDOFMapOutputArray(dof_names, dof_mapping);

    std::vector<std::string> column_names;
    column_names.push_back("x");
    column_names.push_back("y");
    column_names.push_back("z");
    column_names.push_back("id");
    for (auto dof_name_i : dof_names)
        column_names.push_back(dof_name_i);

    std::string base_filename = getFilename("DOF_" + variable_name);
    outputToCSV(dofMapArray, column_names, _num_elems, column_names.size(), base_filename);
}

void
FVDOFMappingOutput::outputAllDOFMaps()
{}

std::map<int, std::map<std::string, Real>>
FVDOFMappingOutput::getSingleVariableDOFMap(std::string variable_name)
{
    std::map<int, std::map<std::string, Real>> dof_map;

    int variable_number = _nl.getVariable(_tid, variable_name).number();
    std::vector<std::string> dof_names;
    std::vector<std::vector<dof_id_type>> dof_mapping;
    getSingleVariableDOFMapping(variable_number, dof_names, dof_mapping);

    Real** dofMapArray;
    dofMapArray = createDOFMapOutputArray(dof_names, dof_mapping);

    std::vector<std::string> column_names;
    column_names.push_back("x");
    column_names.push_back("y");
    column_names.push_back("z");
    column_names.push_back("id");
    for (auto dof_name_i : dof_names)
        column_names.push_back(dof_name_i);
    
    for (int elem_idx : make_range(_num_elems))
    {
        dof_map[elem_idx] = {};
        for (int i=0; i<column_names.size(); i++)
            dof_map[elem_idx][column_names[i]] = dofMapArray[elem_idx][i];
    }
    return dof_map;
}

std::vector<std::string>
FVDOFMappingOutput::getDOFMapColumnNames(std::map<int, std::map<std::string, Real>> dof_map)
{
    std::vector<std::string> column_names;
    for (auto it=dof_map[0].begin(); it!=dof_map[0].end(); it++)
        column_names.push_back(it->first);
    return column_names;
}