#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <MooseTypes.h>
#include "libmesh/petsc_vector.h"

#include "SystemBase.h"
#include "NonlinearSystemBase.h"

#include <fstream>
#include <chrono>

class RomSolverUtilities
{
public:
	// Constructor
	RomSolverUtilities();
	
	// Getters/Setters
	void setSolnAddSolnMult_Vectors(std::string s_add_filepath, std::string s_mult_filepath, const libMesh::Parallel::Communicator * comm);
	void setLogFiles(std::vector<std::string> filenames, bool save, bool timer);
	
	void setNumberNodesRanks(int nNodes_input, int nRanks_input) {nNodes = nNodes_input; nRanks = nRanks_input;}
	
	int getNumNodes() {return nNodes;}
	int getNumRanks() {return nRanks;}
	
	// Setup Functions
	void loadQMatrix(std::string filepath, int nRanks_input, const libMesh::Parallel::Communicator * comm);
	DenseVector<Real> loadVector(std::string filepath, const libMesh::Parallel::Communicator * comm);
	
	// Log-File Functions
	void saveStdMatrixToLog(std::vector<std::vector<Real>> matrix, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter);
	void saveSparseMatrixToLog(SparseMatrix<Real> * matrix, std::string header, std::ofstream & log_file, bool include_iter, bool sparse, int iter);
	void saveDenseMatrixToLog(DenseMatrix<Real> & matrix, std::string header, std::ofstream & log_file, bool include_iter, bool sparse, int iter);
	void saveVectorToLog(DenseVector<Real> vector, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter);
	void saveVectorToLog(const NumericVector<Real> * vector, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter);
	void saveStringToTimeLog(std::string type_name, Real delta_t, std::string units="s") {time_log_file << type_name << ": " << delta_t << " [" << units << "]" << std::endl;};
	void saveStringToTimeLog(std::string type_name, int delta_t, std::string units="s") {time_log_file << type_name << ": " << delta_t << " [" << units << "]" << std::endl;};

	// Timing Functions
	int currentTime(std::string units);
	int currentTime_s() {int now = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count(); return now;};
	int currentTime_ms() {int now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(); return now;};
	int currentTime_us() {int now = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count(); return now;};
	
	// Error Functions
	Real computeRMSE(DenseVector<Real> vec1, DenseVector<Real> vec2);
	Real computeMaxError(DenseVector<Real> vec1, DenseVector<Real> vec2);
	Real computeMaxPercentError(DenseVector<Real> vec1, DenseVector<Real> vec2);
	
	// Mapping Functions
	std::vector<std::vector<Real>> getJacobianRowVectorNonzeroIndices(SparseMatrix<Real> * jacobian);
	std::map<int, std::vector<int>> getJacobianZerosMap(SparseMatrix<Real> * jacobian);
	std::map<int, std::pair<std::string, Real>> getBCMap(SparseMatrix<Number> * jacobian, std::map<std::string, std::vector<dof_id_type>> dof_map, DenseVector<Real> initial_solution);
	std::map<std::string, std::vector<dof_id_type>> getDOFMapping_Old(std::vector<MooseVariableFieldBase * > & variables, MooseMesh & mesh);
	std::map<dof_id_type, std::string> getDOFMapping(std::vector<MooseVariableFieldBase*> & variables, MooseMesh & mesh);
	// Type-Transfer Functions
	DenseMatrix<Number> SparseToDenseMatrix(SparseMatrix<Number> * S);
	DenseMatrix<Number> SparseToDenseMatrix(SparseMatrix<Number> * S, std::map<int, std::vector<int>> matrix_zeros_map);
	DenseVector<Number> NumericToDense(NumericVector<Real> & N);
	
	// Reduction Functions
	DenseMatrix<Number> computeJQ(SparseMatrix<Real> * J, std::map<int, std::vector<int>> row_vector_indicies);
	DenseMatrix<Number> reduceMatrix(SparseMatrix<Real> * M, std::map<int, std::vector<int>> row_vector_indicies);
	DenseVector<Real> reduceVector(DenseVector<Real> & V);
	
	// Solution Update
	DenseVector<Real> reconstructVector(DenseVector<Real> & Vr);
	DenseVector<Real> reconstructOnBoundary(DenseVector<Real> & V, std::map<int, std::pair<std::string, Real>> BC_map);
	DenseVector<Real> updateSolutionOnBoundary(DenseVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map);
	DenseVector<Real> updateSolutionOnBoundary(NumericVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map);
	void updateSolutionOnBoundary_N(NumericVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map);
	
	// Scaling
	DenseVector<Real> scaleVector(DenseVector<Real> V, DenseVector<Real> add_vec, DenseVector<Real> mult_vec, bool inverse, bool do_add=true, bool do_mult=true);
	
protected:
	// Log Files
	bool _save_to_logs;
	bool _save_timing;
	
	std::ofstream gen_log_file;
	std::ofstream sol_log_file;
	std::ofstream res_log_file;
	std::ofstream jac_log_file;
	std::ofstream time_log_file;
	
	DenseMatrix<Real> Q_dense;
	int nRanks;
	int nNodes;
	
	std::map<std::string, std::vector<dof_id_type>> dof_map_old;
	std::map<dof_id_type, std::string> dof_map;

};