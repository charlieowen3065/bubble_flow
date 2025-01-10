#pragma once

#include "FVOutputBase.h"

class FVDOFMappingOutput : public FVOutputBase
{
public:
    static InputParameters validParams();
    FVDOFMappingOutput(const InputParameters & parameters);

    void output();

    std::map<int, std::map<std::string, Real>> getSingleVariableDOFMap(std::string variable_name);
    std::vector<std::string> getDOFMapColumnNames(std::map<int, std::map<std::string, Real>> dof_map);
    
protected:
    dof_id_type** getDOFIndices_SingleVariable(std::vector<dof_id_type> & di, const unsigned int sys_number,
                                               const unsigned int var_number, unsigned int number_componets);
    void getSingleVariableDOFMapping(const unsigned int variable_number, std::vector<std::string> & dof_names, 
                                     std::vector<std::vector<dof_id_type>> & dof_mapping);
    Real** createDOFMapOutputArray(std::vector<std::string> & dof_names,
                                   std::vector<std::vector<dof_id_type>> & dof_mapping);
    void outputSingleDOFMap(std::string variable_name);
    void outputAllDOFMaps();

    
};
