#pragma once

#include "LinearSystemBaseOutput.h"

class DOFMappingOutput : public LinearSystemBaseOutput
{
public:
    static InputParameters validParams();
    DOFMappingOutput(const InputParameters & parameters);

    void output();
protected:
    dof_id_type** getDofIndices_SingleVariable(std::vector<dof_id_type> & di, const unsigned int sys_number,
                                               const unsigned int var_number, unsigned int number_componets);
    void getSingleVariableDOFMapping(const unsigned int variable_number, std::vector<std::string> & dof_names, 
                                     std::vector<std::vector<dof_id_type>> & dof_mapping);
    Real** createDOFMapOutputArray(std::vector<std::string> & dof_names,
                                   std::vector<std::vector<dof_id_type>> & dof_mapping);
    void outputSingleDOFMap(std::string variable_name);
    void outputAllDOFMaps();
};
