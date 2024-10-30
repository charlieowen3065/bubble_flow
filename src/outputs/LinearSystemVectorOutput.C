#include "LinearSystemVectorOutput.h"

registerMooseObjectAliased("BubbleFlowApp", LinearSystemVectorOutput, "LinearSystemVectorOutput");

InputParameters
LinearSystemVectorOutput::validParams()
{
    InputParameters params = LinearSystemBaseOutput::validParams();

    params.addParam<std::string>("vector_type", "Type of vector to output {solution, residual}");

    return params;
}

LinearSystemVectorOutput::LinearSystemVectorOutput(const InputParameters & parameters)
  : LinearSystemBaseOutput(parameters),
    _system_vector(static_cast<ImplicitSystem &>(_nl.system()).get_vector(0)),
    _vector_type(getParam<std::string>("vector_type"))
{
    loadVector();
}

void
LinearSystemVectorOutput::output()  
{
    variableVectors();
}

void
LinearSystemVectorOutput::loadVector()
{
    if (_vector_type == "solution")
        _system_vector = _nl.solution();
    else if (_vector_type == "residual")
        _system_vector = _nl.RHS();
    else
        mooseError("ERROR: " + _vector_type + " is invalid.");
}

void
LinearSystemVectorOutput::outputVectorOfSingleVariable(std::string variable_name)
{
    int variable_number = _nl.getVariable(_tid, variable_name).number();
    unsigned int number_componets = _mesh.nodePtr(0)->n_comp(_sys_number, variable_number);
    std::vector<dof_id_type> dof_indices;
    dof_id_type** component_array = getDOFIndices_SingleVariable(dof_indices, _sys_number,
                                                                 variable_number, number_componets);
    Real** solutionArray = createVectorOutputArrayOfSingleVariable(component_array, number_componets);
    
    // Output
    std::vector<std::string> column_names;
    column_names.push_back("id");
    if (number_componets == 1) {column_names.push_back("u");}
    else {
        std::vector<std::string> vec_components = {"x", "y", "z"};
        for (auto i : make_range(number_componets)) {column_names.push_back("u_" + vec_components[i]);}
    }
    std::string base_filename = getFilename(variable_name + "_" + _vector_type);
    outputToCSV(solutionArray, column_names, 1 + number_componets, _num_nodes, base_filename);
}

Real**
LinearSystemVectorOutput::createVectorOutputArrayOfSingleVariable(dof_id_type** component_array, int number_components)
{
    Real** OutputArray = new Real*[_num_nodes];
    
    for (auto node_idx : make_range(_num_nodes))
    {
        OutputArray[node_idx] = new Real[1 + number_components];
        OutputArray[node_idx][0] = node_idx;
        for (auto ci : make_range(number_components))
        {
            dof_id_type dof_i = component_array[node_idx][4 + ci];
            OutputArray[node_idx][1+ci] = _system_vector.el(dof_i);
        }
    }
    return OutputArray;
}

void
LinearSystemVectorOutput::variableVectors()
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
        outputVectorOfSingleVariable(var_name);
    }
}

void
LinearSystemVectorOutput::fullVector()
{}
