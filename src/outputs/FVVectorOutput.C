#include "FVVectorOutput.h"

registerMooseObjectAliased("BubbleFlowApp", FVVectorOutput, "FVVectorOutput");

InputParameters
FVVectorOutput::validParams()
{
    InputParameters params = FVSystemOutputBase::validParams();

    params.addParam<std::string>("vector_type", "none", "Type of vector to output {none, solution, residual, all}");

    return params;
}

FVVectorOutput::FVVectorOutput(const InputParameters & parameters)
  : FVSystemOutputBase(parameters),
    _vector_type(getParam<std::string>("vector_type"))
{
}

void
FVVectorOutput::output()
{
    // Setup
    setDofMaps();

    // Collects the vector names
    std::vector<std::string> vector_names;
    if (_vector_type == "none") mooseError("Error: please enter a vector name (default == none)");
    else if ((_vector_type == "solution") || (_vector_type == "residual")) vector_names.push_back(_vector_type);
    else if (_vector_type == "all") {vector_names.push_back("solution"); vector_names.push_back("residual");}
    else mooseError("Error: please enter valid vector name ('" + _vector_type + "' is invalid)");

    // Maps vector name to vector reference
    std::map<std::string, NumericVector<Real> *> vector_map;
    vector_map["solution"] = &_nl.solution();
    vector_map["residual"] = &_nl.RHS();

    // Collects the variable names
    std::vector<VariableName> variable_names;
    if ((_output_variable_name != "all") && (_output_variable_name != "none"))
        variable_names.push_back(_output_variable_name);
    else if (_output_variable_name == "all") 
        for (auto v : _nl.getVariableNames())
            variable_names.push_back(v);

    // Outputs vectors
    for (auto var : variable_names)
        for (auto vec : vector_names)
            outputVector(var, vec, "all", *vector_map[vec]);

}

FVVectorTupleMap
FVVectorOutput::extractVariableVector(std::string variable_name, std::string component_name, NumericVector<Real> & vector)
{
    FVVectorTupleMap vector_tuple_map;
    FVDofTupleMap dof_tuple_map = extractDofTuples(variable_name, component_name);
    for (auto it=dof_tuple_map.begin(); it!=dof_tuple_map.end(); it++){
        vector_tuple_map[it->first] = 
            std::make_tuple(std::get<0>(it->second),
                            std::get<1>(it->second),
                            vector.el(it->first));
    }
    return vector_tuple_map;
}

void
FVVectorOutput::outputVector(std::string variable_name, std::string vector_name, 
                             std::string component_name, NumericVector<Real> & vector)
{
    FVVectorTupleMap vector_tuple_map = extractVariableVector(variable_name, component_name, vector);
    vectorMapToCSV(vector_tuple_map, vector_name + "_" + variable_name);
}