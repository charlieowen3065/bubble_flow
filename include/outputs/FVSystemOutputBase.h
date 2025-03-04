#pragma once

#include "FVOutputBase.h"
#include "FVDOFMappingOutput.h"

class FVSystemOutputBase : public FVOutputBase
{
public:
    static InputParameters validParams();
    FVSystemOutputBase(const InputParameters & parameters);

    void output();
protected:
    // -------- Class Functions -------- //
    // General setup functions
    void setupTimeVariables();
    void setVariablesList();
    void setDofMaps();
    bool doOutput();

    // Extraction functions
    FVDofTupleMap extractDofTuples(std::string variable_name, std::string component_name);
    std::pair<std::vector<std::string>, std::vector<std::string>> extractPossibleVectorComponents(std::vector<std::string> column_names, std::string variable_name);

    // -------- Class Variables -------- //
    FVDOFMappingOutput dof_object;

    // Timing Inputs
    std::string _which_time;
	std::string _time_units;
    std::vector<int> _timesteps_to_display;

    /// Latter Initilizations
    // setupTimeVariables
    bool _show_converged_time;
	bool _show_nl_time;
	bool _show_l_time;

    // Variable list
    std::vector<std::string> varaible_names;

    // Dof Maps
    std::map<std::string, FVDofMap> dof_maps;
    std::map<std::string, std::vector<std::string>> dof_column_names;

    // START-TIME
    int START_TIME_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

};
