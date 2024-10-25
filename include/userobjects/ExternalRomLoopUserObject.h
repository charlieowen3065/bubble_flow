#pragma once

#include "RomSolverUtilities.h"
#include "GeneralUserObject.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_linear_solver.h"

class ExternalRomLoopUserObject : public GeneralUserObject, public RomSolverUtilities
{
public:
	// Constructor, etc.
	static InputParameters validParams();
	ExternalRomLoopUserObject(const InputParameters & parameters);
	virtual ~ExternalRomLoopUserObject();
	
	// LOOPS
	void FOMLoop();
	void FOMLoopWithLinear();
	void FOMLoopDEBUG();
	void ROMLoop();
	void ROMLoopDEBUG();
	
	// Solution Update
	void updateSolution(DenseVector<Real> new_solution);
	void updateSolution(NumericVector<Number> * new_solution);
	
	// Line Search Algorithms
	DenseVector<Real> lineSearchAlgorithm1(DenseMatrix<Real> & J);
	Real getOptimalAlpha(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0);
	DenseVector<Real> LinearUpdateFOM(DenseMatrix<Real> * J);
	DenseVector<Real> LinearUpdateROM(DenseMatrix<Real> * Jr);
	
	Real changeAlpha_Type1(Real alpha, Real beta);
	Real changeAlpha_Type2(Real alpha, Real beta, Real error_ratio, int iter);
	Real changeAlpha_Type3(Real alpha, Real beta, Real err, Real err0, Real err_prev, Real err_prev2, int iter);
	
	Real alphaComputeGradient(Real alpha, PetscVector<Number> & soln_inp, PetscVector<Number> & del_soln);
	Real alphaOptimizationGradientDescent(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0);
	std::pair<Real, std::vector<Real>> alphaOptimizationGradientDescent2(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0);
	DenseVector<Real> alphaOptimizationGradientDescent3(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0);
	void alphaOptimizationGradientDescent4(NumericVector<Real> * curr_soln, DenseVector<Real> u, DenseVector<Real> delu, Real alpha0);
	
	
	Real getStepSize(std::vector<Real> error_history, Real input_step_size);
	
	// Startup
	void initialSetup_ROM();
	void initializeDOFMapping();
	void initializeLogFiles();
	void initializeQMatrix();
	void initializeConvergenceSettings();
	
	void initalizeSolutionScaling();
	void initalizeSolutionScalingVectors();
	
	// Initalization
	void initializeMatrixAndVectors();
	void initializeReferenceSolution();
	void initializeBoundaryConditions();
	void initializeConstantBoundaryConditions();
	void initializeFunctionBoundaryConditions();
	void initializeCSVBoundaryConditions();

protected:
	// Required functions for user-objects
	virtual void initialize();
	virtual void finalize();
	virtual void execute();
	
	bool _perform_startup_initalization = true;
	
	/// Matrix and Vector initializations
	SparseMatrix<Number> * _jac_matrix;
	NumericVector<Number> * _residual;
	NumericVector<Number> * _curr_sol;
	
	/// User Inputs
	// System Parameters
	NonlinearSystemBase & _nl_sys;
	MooseMesh & _mesh;
	THREAD_ID _tid;
	// Control Booleans
	bool _save_to_logs;
	bool _save_timing;
	bool _save_all_errors;
	bool _save_final_solution;
	bool _use_rom_loop;
	bool _transient;
	bool _fom_debug;
	bool _optimize_alpha_basic;
	bool _optimize_alpha_grad_descent;
	bool _use_linear_iterations;
	// Convergence Inputs
	int _convergence_type;
	int _alpha_optimization_type;
	Real _abs_tol_rom;
	Real _rmse_tol_rom;
	Real _max_err_tol_rom;
	Real _max_per_err_tol_rom;
	int _max_iter_rom;
	
	Real _alpha_tol_ROM;
	int _max_alpha_iter;
	Real _alpha_opt_beta;
	Real _linear_tol_ROM;
	int _max_linear_iter;
	
	// CSV Paths
	std::string _Qpath_rom;
	std::string _expected_soln_path;
	// Scaling Parameters
	Real _relaxation_factor;
	bool _scaling_vectors;
	std::string _addition_vector_csv_path;
	std::string _multiplication_vector_csv_path;
	// Boundary Conditions
	bool _function_BC;
	bool _csv_BC;
	std::string _BC_filepath;
	
	/// Log Files
	std::ofstream gen_log_file;        // General Outputs
	std::ofstream jac_log_file;        // Jacobian Data
	std::ofstream res_log_file;        // Residual Data
	std::ofstream sol_log_file;        // Solution Data
	std::ofstream err_log_file;        // Refernce Errors
	std::ofstream all_err_log_file;    // All Error Types
	std::ofstream conv_sol_log_file;   // Converged Solutions
	std::ofstream final_sol_log_file;  // All Final NL Solutions
	std::ofstream time_log_file;       // Timings
	
	std::ofstream debug_log_file;
	
	/// Initalization of Convergence Parameters
	int iter;  // Iteration number
	Real _tolerance;  // Overwritten with desired tolerance
	std::string _convergence_name;
	
	Real _dt_rom;
	
	/// Initalization of Q-Matrix Related Attributes
	int nRanks;
	int _nRanks_inp;
	int nNodes;
	
	/// Initalization of Expected-Solution Related Attributes
	DenseVector<Real> _expected_soln;
	bool _scale_solution;
	DenseVector<Real> _addition_scaling_vector;
	DenseVector<Real> _multiplication_scaling_vector;
	
	/// Mappings
	std::map<int, std::pair<std::string, Real>> BC_map;
	std::vector<std::vector<Real>> row_vector_indicies;
	std::map<int, std::vector<int>> jacobian_zeros_map;
	std::map<dof_id_type, std::string> dof_map;
	int num_dofs;
	
	
	// Linear Solver
	PetscLinearSolver<Real> linear_solver;
	
};