#pragma once

#include "FileOutput.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"
#include "FEProblemBase.h"

#include <chrono>
#include <ctime>

// Type definitions
typedef std::map<int/*node-idx*/, std::map<std::string/*column-name*/, Real/*dof-#*/>> FVDofMap;
typedef std::map<int/*dof-#*/, std::tuple<int/*node-idx*/, std::string/*vec-comp*/>> FVDofTupleMap;
typedef std::map<int/*dof-#*/, std::tuple<int/*node-idx*/, std::string/*vec-comp*/, Real/*value*/>> FVVectorTupleMap;

class FVOutputBase : public FileOutput
{
public:
    static InputParameters validParams();
    FVOutputBase(const InputParameters & parameters);

    // Main Output
    void output();

protected:
    // -------- Class Functions -------- //
    // Output functions
    std::string getFilename(std::string BASE);
    void outputToCSV(Real** arr_2d, std::vector<std::string> column_names, int num_rows, int num_cols, std::string base_name);
    void vectorMapToCSV(FVVectorTupleMap vector_map, std::string base_filename);

    // Spatial functions
    std::vector<Real> getNodeCoordinates(const Point & p, const Real & id);
    std::vector<Real> getElementCoordinates(const Elem & elem, const Real & id);

    // Boolean functions
    bool isInVector(std::vector<std::string> vec, std::string target);
    
    // -------- Class Variables -------- //
    /// Class Initilization
    // Nonlinear system variables
    THREAD_ID _tid;
    const unsigned int _nl_sys_num;
	NonlinearSystemBase & _nl;
    
    // Nonlinear-Variable Information
	std::string _output_variable_name;
    
    // System variables
    const System & _sys;
	const DofMap & _dof_map;
	const unsigned int _sys_number;
	std::string _system_name;
    
    // Meshing
    MooseMesh & _mesh;
	dof_id_type _num_nodes;
    dof_id_type _num_elems;
    
};