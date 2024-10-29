#include "SolutionVectorOutput.h"

registerMooseObjectAliased("BubbleFlowApp", SolutionVectorOutput, "SolutionVectorOutput");

InputParameters
SolutionVectorOutput::validParams()
{
    InputParameters params = LinearSystemBaseOutput::validParams();

    return params;
}

SolutionVectorOutput::SolutionVectorOutput(const InputParameters & parameters)
  :  LinearSystemBaseOutput(parameters)
{}

void
SolutionVectorOutput::output()  
{
    variableSolutions();
}

void
SolutionVectorOutput::outputSolutionOfSingleVariable(std::string variable_name)
{
    int variable_number = _nl.getVariable(_tid, variable_name).number();
    unsigned int number_componets = _mesh.nodePtr(0)->n_comp(_sys_number, variable_number);
    std::vector<dof_id_type> dof_indices;
    dof_id_type** component_array = getDOFIndices_SingleVariable(dof_indices, _sys_number,
                                                                 variable_number, number_componets);
    Real** solutionArray = createSolutionOutputArrayOfSingleVariable(component_array, number_componets);
    
    // Output
    std::vector<std::string> column_names;
    column_names.push_back("id");
    if (number_componets == 1) {column_names.push_back("u");}
    else {
        std::vector<std::string> vec_components = {"x", "y", "z"};
        for (auto i : make_range(number_componets)) {column_names.push_back("u_" + vec_components[i]);}
    }
    std::string base_filename = getFilename(variable_name + "_solution");
    outputToCSV(solutionArray, column_names, 1 + number_componets, _num_nodes, base_filename);
}

Real**
SolutionVectorOutput::createSolutionOutputArrayOfSingleVariable(dof_id_type** component_array, int number_components)
{
    Real** OutputArray = new Real*[_num_nodes];
    const NumericVector<Real> * current_solution = _nl.currentSolution();
    
    for (auto node_idx : make_range(_num_nodes))
    {
        OutputArray[node_idx] = new Real[1 + number_components];
        OutputArray[node_idx][0] = node_idx;
        for (auto ci : make_range(number_components))
        {
            dof_id_type dof_i = component_array[node_idx][4 + ci];
            OutputArray[node_idx][1+ci] = current_solution->el(dof_i);
        }
    }
    return OutputArray;
}

void
SolutionVectorOutput::variableSolutions()
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
        outputSolutionOfSingleVariable(var_name);
    }
}

void
SolutionVectorOutput::fullSolution()
{}
