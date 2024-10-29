#include "DOFMappingOutput.h"

registerMooseObjectAliased("BubbleFlowApp", DOFMappingOutput, "DOFMappingOutput");

InputParameters
DOFMappingOutput::validParams()
{
    InputParameters params = LinearSystemBaseOutput::validParams();

    return params;
}

DOFMappingOutput::DOFMappingOutput(const InputParameters & parameters)
  : LinearSystemBaseOutput(parameters)
{}

void
DOFMappingOutput::output() 
{
    std::vector<VariableName> var_names;
    if ((_output_variable_name != "all") && (_output_variable_name != "none"))
    {
        var_names.push_back(_output_variable_name);
    } else if (_output_variable_name == "all")
    {
        for (auto v : _nl.getVariableNames())
            var_names.push_back(v);
    }

    for (auto var_name : var_names){
        outputSingleDOFMap(var_name);
    }
}

void
DOFMappingOutput::getSingleVariableDOFMapping(const unsigned int variable_number,
                                              std::vector<std::string> & dof_names, 
                                              std::vector<std::vector<dof_id_type>> & dof_mapping)
{
	bool is_vector;
    unsigned int number_componets = _mesh.nodePtr(0)->n_comp(_sys_number, variable_number);
    if (number_componets == 1) {is_vector = false;} else {is_vector = true;}

	std::vector<dof_id_type> dof_indices_full;
	dof_id_type** component_array = getDOFIndices_SingleVariable(dof_indices_full, _sys_number,
                                                                 variable_number, number_componets);
    
    if (!is_vector)
    {
        dof_names.push_back("dof");
        std::vector<dof_id_type> dof_map_i;
        for (auto i : dof_indices_full)
            dof_map_i.push_back(i);
        dof_mapping.push_back(dof_map_i);
    } 
    else 
    {
        std::vector<std::string> vec_comps = {"x", "y", "z"};
        for (auto c_i : make_range(number_componets))
		{
            dof_names.push_back("dof_" + vec_comps[c_i]);
			std::vector<dof_id_type> dof_map_i;
			for (auto i : make_range(_num_nodes)) 
				dof_map_i.push_back(component_array[c_i][i]);
			dof_mapping.push_back(dof_map_i);		
		}
    }
}

Real**
DOFMappingOutput::createDOFMapOutputArray(std::vector<std::string> & dof_names,
                                          std::vector<std::vector<dof_id_type>> & dof_mapping)
{
    Real** OutputArray = new Real*[_num_nodes];
    int _number_componets = dof_names.size();

    for (dof_id_type node_idx : make_range(_num_nodes))
    {
        const Node * node_i = _mesh.nodePtr(node_idx);
        std::vector<Real> coords = getNodeCoordinates(*node_i, node_idx);

        OutputArray[node_idx] = new Real[4 + _number_componets];
        OutputArray[node_idx][0] = coords[0];
        OutputArray[node_idx][1] = coords[1];
        OutputArray[node_idx][2] = coords[2];
        OutputArray[node_idx][3] = coords[3];
        
        for (auto ci : make_range(_number_componets))
        {
            dof_id_type dof_i = dof_mapping[ci][node_idx];
            OutputArray[node_idx][4 + ci] = dof_i;
        }
    }
    return OutputArray;
}

void
DOFMappingOutput::outputSingleDOFMap(std::string variable_name)
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

    std::string base_filename = getFilename(variable_name);
    outputToCSV(dofMapArray, column_names, _num_nodes, column_names.size(), base_filename);
}
