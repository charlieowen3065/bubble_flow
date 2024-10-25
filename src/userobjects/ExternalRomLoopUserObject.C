#include "ExternalRomLoopUserObject.h"
#include "DelimitedFileReader.h"

registerMooseObject("MalamuteApp", ExternalRomLoopUserObject);

InputParameters
ExternalRomLoopUserObject::validParams()
{
	InputParameters params = GeneralUserObject::validParams();
	
	/// Control Booleans
	// Log Files
	params.addParam<bool>("save_to_logs", true, "Save Jacobian, Solution, and Residual data to log files");
	params.addParam<bool>("save_timing_log", true, "Save timing data to corresponding log file");
	params.addParam<bool>("save_all_error_types", true, "Saves all error types, not just the reference error");
	params.addParam<bool>("save_final_solution", true, "Saves all error types, not just the reference error");
	// ROM Specific
	params.addParam<bool>("use_rom_loop", false, "Boolean to control if the ROM loop is used, or the FOM loop");
	// Run Type
	params.addParam<bool>("transient", false, "Boolean to control if the loop is transient or steady-state");
	params.addParam<bool>("fom_debug", false, "Boolean to control the FOM debug loop");
	params.addParam<bool>("optimize_alpha_basic", false, "Boolean to control if 'alpha' is opitimized");
	params.addParam<bool>("optimize_alpha_grad_descent", false, "Boolean to control if 'alpha' is opitimized");
	params.addParam<bool>("use_linear_iterations", false, "Boolean to control if linear iterations are used");
	
	/// Convergance Inputs
	// Type control
	params.addParam<int>("convergence_type", 0, "Type of convergence {0, -1 : Absolute Tolerance, \
																	  1     : Root-Mean Squared Error, \
																	  2     : Maximum Error, \
																	  3     : Percent Error}");
	params.addParam<int>("alpha_optimization_type", 0, "Alpha optimization algorithms {0, 1}");
	// Tolerance
	params.addParam<Real>("abs_tol_ROM", 1.0, "Absolute tolerance (residual-based)");
	params.addParam<Real>("rmse_tol_ROM", 1.0, "RMSE tolerance");
	params.addParam<Real>("max_err_tol_ROM", 1.0, "Max Error tolerance");
	params.addParam<Real>("max_per_err_tol_ROM", 1.0, "Percent Error tolerance");
	// Iterations
	params.addParam<int>("max_iter_rom", 50, "max iterations");	
	// Alpha Optimization
	params.addParam<Real>("alpha_tol_ROM", 1.0, "Alpha-optimization tolerance");
	params.addParam<int>("max_alpha_iter", 50, "max iterations for alpha-optimization");
	params.addParam<Real>("alpha_opt_beta", 0.5, "Alpha-optimization beta value");	
	// Linear Iterations
	params.addParam<Real>("linear_tol_ROM", 1.0, "Alpha-optimization tolerance");
	params.addParam<int>("max_linear_iter", 50, "max iterations for alpha-optimization");
	
	/// CSV Paths
	params.addParam<std::string>("Q_path", "", "Path to Q-Matrix");
	params.addParam<std::string>("expected_solution_path", "", "Path to the expected solution");
	params.addParam<int>("nRanks_inp", -1, "number of ranks desired (-1 == all ranks)");
	
	/// Scaling Parameters
	params.addParam<Real>("relaxation_factor", 1.0, "Relaxation factor for .add() function");
	
	params.addParam<bool>("scaling_vectors", false, "Boolean for calling the scaling vector csv files");
	params.addParam<std::string>("addition_vector_csv_path", "", "csv filepath for the addition vector");
	params.addParam<std::string>("multiplication_vector_csv_path", "", "csv filepath for the multiplication vector");
	
	/// Boundary Conditions
	params.addParam<bool>("function_BC", false, "Inputs wether the BCs are some function, or are constant through time/states");
	params.addParam<bool>("csv_BC", false, "Inputs wether the BCs are some CSV, or are constant through time/states");
	params.addParam<std::string>("BC_filepath", "Filepath to the BC CSV file");
	
	return params;
}

ExternalRomLoopUserObject::ExternalRomLoopUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
	RomSolverUtilities(),
	/// System Parameters
	_nl_sys(_fe_problem.getNonlinearSystemBase(0)),
	_mesh(_fe_problem.mesh()),
	_tid(getParam<THREAD_ID>("_tid")),
	/// Control Booleans
	_save_to_logs(getParam<bool>("save_to_logs")),
	_save_timing(getParam<bool>("save_timing_log")),
	_save_all_errors(getParam<bool>("save_all_error_types")),
	_save_final_solution(getParam<bool>("save_final_solution")),
	_use_rom_loop(getParam<bool>("use_rom_loop")),
	_transient(getParam<bool>("transient")),
	_fom_debug(getParam<bool>("fom_debug")),
	_optimize_alpha_basic(getParam<bool>("optimize_alpha_basic")),
	_optimize_alpha_grad_descent(getParam<bool>("optimize_alpha_grad_descent")),
	_use_linear_iterations(getParam<bool>("use_linear_iterations")),
	/// Convergance Inputs
	_convergence_type(getParam<int>("convergence_type")),
	_alpha_optimization_type(getParam<int>("alpha_optimization_type")),
	_abs_tol_rom(getParam<Real>("abs_tol_ROM")),
	_rmse_tol_rom(getParam<Real>("rmse_tol_ROM")),
	_max_err_tol_rom(getParam<Real>("max_err_tol_ROM")),
	_max_per_err_tol_rom(getParam<Real>("max_per_err_tol_ROM")),
	_max_iter_rom(getParam<int>("max_iter_rom")),
	
	_alpha_tol_ROM(getParam<Real>("alpha_tol_ROM")),
	_max_alpha_iter(getParam<int>("max_alpha_iter")),
	_alpha_opt_beta(getParam<Real>("alpha_opt_beta")),
	_linear_tol_ROM(getParam<Real>("linear_tol_ROM")),
	_max_linear_iter(getParam<int>("max_linear_iter")),
	/// CSV Paths
	_Qpath_rom(getParam<std::string>("Q_path")),
	_expected_soln_path(getParam<std::string>("expected_solution_path")),
	/// Scaling Parameters
	_relaxation_factor(getParam<Real>("relaxation_factor")),
	_scaling_vectors(getParam<bool>("scaling_vectors")),
	_addition_vector_csv_path(getParam<std::string>("addition_vector_csv_path")),
	_multiplication_vector_csv_path(getParam<std::string>("multiplication_vector_csv_path")),
	/// Boundary Conditions
	_function_BC(getParam<bool>("function_BC")),
	_csv_BC(getParam<bool>("csv_BC")),
	_BC_filepath(getParam<std::string>("BC_filepath")),
	/// Log File Initalization
	gen_log_file("general_log_file.txt"),
	jac_log_file("jacobian_log_file.txt"),
	res_log_file("residual_log_file.txt"),
	sol_log_file("solution_log_file.txt"),
	err_log_file("error_log_file.txt"),
	all_err_log_file("all_errors_log_file.txt"),
	conv_sol_log_file("converged_solution_log_file.txt"),
	final_sol_log_file("final_solution_log_file.txt"),
	time_log_file("timing_log_file.txt"),
	debug_log_file("debug_log_file.txt"),
	// Convergence Parameters
	iter(0),
	// Q-based input
	_nRanks_inp(getParam<int>("nRanks_inp")),
	// Solver
	linear_solver(_communicator)
{
	// Checks
	if ((_csv_BC) and (_BC_filepath.empty()))
		mooseError("Please input B.C. filepath and corresponding boolean vector");
	
	
}

ExternalRomLoopUserObject::~ExternalRomLoopUserObject() {}

// ************************* Main Functions ************************* //
void ExternalRomLoopUserObject::initialize() 
{
	std::cout << "***** HERE INITALIZE *****" << std::endl;
	std::cout << "_t_step: " << _t_step << std::endl;
	std::cout << "_t: " << _t << std::endl;
	std::cout << "_dt: " << _dt << std::endl;
	
	initializeMatrixAndVectors();
	if (_perform_startup_initalization)
		initialSetup_ROM();
	initializeBoundaryConditions();
	if (_convergence_type >= 0)
		initializeReferenceSolution();
}
void ExternalRomLoopUserObject::finalize() 
{
	std::cout << "***** HERE FINAL *****" << std::endl;
	std::cout << "_t_step: " << _t_step << std::endl;
	std::cout << "_t: " << _t << std::endl;
	std::cout << "_dt: " << _dt << std::endl;
}

void
ExternalRomLoopUserObject::execute()
{
	if (_use_rom_loop)
		if (!_fom_debug)
			ROMLoop();
		else
			ROMLoopDEBUG();
	else
		if (!_fom_debug)
			FOMLoop();
		else
			FOMLoopDEBUG();
}

// ************************* Startup Functions ************************* //

void
ExternalRomLoopUserObject::initialSetup_ROM()
{
	initializeDOFMapping();
	initializeLogFiles();
	if (_use_rom_loop)
		initializeQMatrix();
	else{
		nNodes = _curr_sol->size();
		nRanks = 1;
		setNumberNodesRanks(nNodes, nRanks);
		}
	initializeConvergenceSettings();
	initalizeSolutionScaling();
	_perform_startup_initalization = false;
}

void
ExternalRomLoopUserObject::initializeDOFMapping()
{
	std::vector<MooseVariableFieldBase*> variables = _nl_sys.getVariables(_tid);
	dof_map = getDOFMapping(variables, _mesh);
	num_dofs = dof_map.size();
}

void
ExternalRomLoopUserObject::initializeLogFiles()
{
	std::vector<std::string> log_files;
	log_files.push_back("general_log_file.txt");
	log_files.push_back("solution_log_file.txt");
	log_files.push_back("residual_log_file.txt");
	log_files.push_back("jacobian_log_file.txt");
	log_files.push_back("timing_log_file.txt");
	setLogFiles(log_files, _save_to_logs, _save_timing);
}

void
ExternalRomLoopUserObject::initializeQMatrix()
{
	loadQMatrix(_Qpath_rom, _nRanks_inp, &_communicator);
	nNodes = getNumNodes();
	nRanks = getNumRanks();
}

void
ExternalRomLoopUserObject::initializeConvergenceSettings()
{
	if (_convergence_type <= 0)
		_tolerance = _abs_tol_rom;
	else if (_convergence_type == 1)
	{
		_tolerance = _rmse_tol_rom;
		_convergence_name = "RMSE";
	}
	else if (_convergence_type == 2)
	{
		_tolerance = _max_err_tol_rom;
		_convergence_name = "MaxError";
	}
	else if (_convergence_type == 3)
	{
		_tolerance = _max_per_err_tol_rom;
		_convergence_name = "MaxPercentError";
	}
	else
		mooseError("Error: Please enter proper convergence_type {-1, 0, 1, 2, 3}");
}

void
ExternalRomLoopUserObject::initalizeSolutionScaling()
{
	if (_scaling_vectors)
		initalizeSolutionScalingVectors();
	
	if (_scaling_vectors)
		_scale_solution = true;
	else
		_scale_solution = false;
}

void
ExternalRomLoopUserObject::initalizeSolutionScalingVectors()
{
	_addition_scaling_vector.resize(nNodes);
	if (_addition_vector_csv_path != "")
		_addition_scaling_vector = loadVector(_addition_vector_csv_path, &_communicator);
	else
		_addition_scaling_vector.zero();
	
	_multiplication_scaling_vector.resize(nNodes);
	if (_multiplication_vector_csv_path != "")
		_multiplication_scaling_vector = loadVector(_multiplication_vector_csv_path, &_communicator);
	else
	{
		_multiplication_scaling_vector.zero();
		for (int i=0; i<nNodes; i++)
			_multiplication_scaling_vector(i) = 1.0;
	}
}

// ************************* Initalization Functions ************************* //
void 
ExternalRomLoopUserObject::initializeMatrixAndVectors()
{
	// Residual & Solution
	_residual = &_nl_sys.RHS();
	_curr_sol = &_nl_sys.solution();
	
	// Assemble matrix & Petsc stuff
	auto & petsc_options = _fe_problem.getPetscOptions();
	auto & pars = _fe_problem.solverParams();
	Moose::PetscSupport::petscSetOptions(petsc_options, pars);
	
	_jac_matrix = &static_cast<ImplicitSystem &>(_nl_sys.system()).get_system_matrix();

	jacobian_zeros_map = getJacobianZerosMap(_jac_matrix);
}

void
ExternalRomLoopUserObject::initializeReferenceSolution()
{
	MooseUtils::DelimitedFileReader soln_csv_reader(_expected_soln_path, &_communicator);
	soln_csv_reader.read();
	const std::vector<std::vector<double>> & soln_inp/*[col][row]*/ = soln_csv_reader.getData();
	
	/*
	_expected_soln.resize(num_dofs);
	
	if (!_transient)
		for (int i=0; i<num_dofs; i++)
			_expected_soln(i) = soln_inp[0][i];
	else
		for (int i=0; i<num_dofs; i++)
			_expected_soln(i) = soln_inp[_t_step-1][i];
	*/
	_expected_soln.resize(nNodes);
	
	if (!_transient)
		for (int i=0; i<nNodes; i++)
			_expected_soln(i) = soln_inp[0][i];
	else
		for (int i=0; i<nNodes; i++)
			_expected_soln(i) = soln_inp[_t_step-1][i];
	
}

void
ExternalRomLoopUserObject::initializeBoundaryConditions()
{
	if (_function_BC)
		initializeFunctionBoundaryConditions();
	if (_csv_BC)
		initializeCSVBoundaryConditions();
	else
		initializeConstantBoundaryConditions();	
}

void
ExternalRomLoopUserObject::initializeConstantBoundaryConditions()
{
	mooseError("Error: Have not implemented ability to constant BCs");
}

void
ExternalRomLoopUserObject::initializeFunctionBoundaryConditions()
{
	mooseError("Error: Have not implemented ability to handle function BCs");
}

void
ExternalRomLoopUserObject::initializeCSVBoundaryConditions()
{
	MooseUtils::DelimitedFileReader BC_csv_reader(_BC_filepath, &_communicator);
	BC_csv_reader.read(); 
	const std::vector<std::vector<double>> & BC_inp/*[col][row]*/ = BC_csv_reader.getData();
	
	for (int i=0; i<num_dofs; i++)
	{
		if (BC_inp[0][i] == 1)
		{
			std::string v_name = dof_map[i];
			Real BC_val = BC_inp[1][i];
			BC_map[i] = std::make_pair(v_name, BC_val);
		}
	}
}

// ************************* Solution Update ************************* //
void
ExternalRomLoopUserObject::updateSolution(DenseVector<Number> new_solution)
{
	DenseVector<Number> S = updateSolutionOnBoundary(new_solution, BC_map);
	NumericVector<Number> * local_soln = &_nl_sys.solution();
	
	for (auto i=0; i<nNodes; i++)
	{
		local_soln->set(i, S.el(i));
	}
	local_soln->close();
	_nl_sys.update();

}

void
ExternalRomLoopUserObject::updateSolution(NumericVector<Number> * new_solution)
{
	DenseVector<Number> S = updateSolutionOnBoundary(*new_solution, BC_map);
	NumericVector<Number> * local_soln = &_nl_sys.solution();
	
	for (auto i=0; i<nNodes; i++)
	{
		local_soln->set(i, S.el(i));
	}
	local_soln->close();
	_nl_sys.update();

}

// ************************* Line Search Algorithms ************************* //
Real
ExternalRomLoopUserObject::getOptimalAlpha(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0)
{
	PetscVector<Number> soln(_communicator, u.size());
	PetscVector<Number> del_soln(_communicator, u.size());
	PetscVector<Number> res(_communicator, u.size());
	for (int i=0; i<u.size(); i++) {
		soln.set(i, u.el(i));
		del_soln.set(i, delu.el(i));
		res.set(i, 0.0);}
	
	Real tol = _alpha_tol_ROM;
	
	_fe_problem.computeResidual(soln, res, 0);
	Real err0 = res.l2_norm();
	
	Real err;
	Real err_prev = std::numeric_limits<double>::max() / 10000.0;
	Real err_prev2 = std::numeric_limits<double>::max() / 10000.0;
	Real alpha = alpha0;
	Real beta = _alpha_opt_beta;
	int max_iter = _max_alpha_iter;
	int iter = 1;
	
	Real min_err_alpha = alpha0;
	Real min_err = std::numeric_limits<double>::max();
	
	bool is_converged = false;
	int const_counter = 0;
	
	std::cout << std::endl;
	
	while ((!is_converged) and (iter < max_iter+1))
	{
		PetscVector<Number> soln_i(_communicator, u.size());
		soln_i.add(1, soln);
		soln_i.add(alpha, del_soln);
		_fe_problem.computeResidual(soln_i, res, 0);
		err = res.l2_norm();
		res.zero();
		soln_i.zero();
		if (err < min_err){
			min_err = err;
			min_err_alpha = alpha;
			const_counter = 0;}

		if (err - err_prev < min_err / 10.0)
			const_counter++;
		else
			const_counter = 0;
		
		if ((err < tol) or (5.0*min_err < err) or (const_counter==10))
			is_converged = true;
		
		std::cout << "        " << iter << " | Alpha: " << alpha << " | Error: " << err;
		std::cout << "  [ tol: " << tol << " | min_err: " << min_err << " | 5.0*min_err: " << 5.0*min_err << " | const_counter: " << const_counter << " ]" << std::endl;
		
		if (_alpha_optimization_type == 0)
			alpha = changeAlpha_Type1(alpha, beta);
		else if (_alpha_optimization_type == 1)
			alpha = changeAlpha_Type2(alpha, beta, err / err0, iter);
		else
			mooseError("ERROR: Please enter valid 'alpha_optimization_type'");
		
		if (iter > 1)
			err_prev2 = err_prev;
		err_prev = err;
		iter++;
	}
	
	return min_err_alpha;
}

Real
ExternalRomLoopUserObject::changeAlpha_Type1(Real alpha, Real beta)
{
	return alpha * beta;
}

Real
ExternalRomLoopUserObject::changeAlpha_Type2(Real alpha, Real beta, Real error_ratio, int iter)
{
	if (error_ratio > 1)
		return alpha / (beta * iter);
	else
		return alpha * (beta * iter);
}

Real
ExternalRomLoopUserObject::alphaOptimizationGradientDescent(DenseVector<Real> u, DenseVector<Real> delu, Real alpha0)
{
	PetscVector<Number> soln(_communicator, u.size());
	PetscVector<Number> del_soln(_communicator, u.size());
	PetscVector<Number> res(_communicator, u.size());
	for (int i=0; i<u.size(); i++) {
		soln.set(i, u.el(i));
		del_soln.set(i, delu.el(i));
		res.set(i, 0.0);}
	
	bool is_converged = false;
	Real tol = _alpha_tol_ROM;
	int max_iter = _max_alpha_iter;
	int iter = 1;
	Real alpha = alpha0;
	
	Real min_err_alpha = alpha0;
	Real min_err = std::numeric_limits<double>::max();
	
	Real grad_f;
	Real step_size = _alpha_opt_beta;
	
	PetscVector<Number> soln_i(_communicator, u.size());
	soln_i.add(1, soln);
	soln_i.add(alpha, del_soln);
	updateSolutionOnBoundary_N(soln_i, BC_map);
	_fe_problem.computeResidual(soln_i, res, 0);
	Real err = res.l2_norm();
	
	if (err < min_err){min_err = err; min_err_alpha = alpha;}
	
	std::cout << "    [Alpha-Opt.] Iteration " << iter << " (Alpha = " << alpha << ")" << std::endl;
	std::cout << "        Residual: " << err << std::endl;
	err_log_file << "    [Alpha-Opt.] Iteration " << iter << " (Alpha = " << alpha << ")" << std::endl;
	err_log_file << "        Residual: " << err << std::endl;
	
	res.zero();
	soln_i.zero();
	
	Real base_grad_f = std::abs(alphaComputeGradient(alpha, soln_i, del_soln)) * 10.0;
	if (base_grad_f == 0)
		base_grad_f = 1.0;
	
	//std::cout << std::endl << "        " << iter-1 << " | Alpha: " << alpha << " | Error: " << err << " | Base-Grad: " << base_grad_f << std::endl;
	
	std::vector<Real> grad_history;
	grad_history.push_back(base_grad_f);
	Real del_alpha = step_size;
	
	while ((!is_converged) and (iter < max_iter+1))
	{
		PetscVector<Number> soln_i(_communicator, u.size());
		soln_i.add(1, soln);
		grad_f = alphaComputeGradient(alpha, soln_i, del_soln);
		del_alpha = getStepSize(grad_history, del_alpha);
		alpha = alpha - del_alpha*(grad_f/std::abs(grad_f));
		
		soln_i.add(alpha, del_soln);
		updateSolutionOnBoundary(soln_i, BC_map);
		//std::cout << "        (Post-" << iter << ") soln_i: [2] = " << soln_i.el(2) << " | [16] = " << soln_i.el(16) << " | [3] = " << soln_i.el(3) << " | [13210] = " << soln_i.el(13210) << std::endl;
		_fe_problem.computeResidual(soln_i, res, 0);
		err = res.l2_norm();
		res.zero();
		soln_i.zero();
		
		grad_history.push_back(grad_f);
		
		if (err < min_err){
			min_err = err;
			min_err_alpha = alpha;}
		
		//std::cout << "        " << iter << " | Alpha: " << alpha << " | Error: " << err << " | del_alpha: " << del_alpha;
		//std::cout << "  [ tol: " << tol << " | grad_f: " << grad_f << " | min_err: " << min_err << " | min_err_alpha: " << min_err_alpha << "]" << std::endl;
		
		bool condition_1 = err < tol;
		bool condition_2 = std::abs(grad_f) < tol;
		bool condition_3 = grad_f == 0;
		bool condition_4 = std::abs(del_alpha / alpha) < 0.0001;
		
		std::cout << "    [Alpha-Opt.] Iteration " << iter << " (Alpha = " << alpha << ")" << std::endl;
		std::cout << "        Residual: " << err << std::endl;
		err_log_file << "    [Alpha-Opt.] Iteration " << iter << " (Alpha = " << alpha << ")" << std::endl;
		err_log_file << "        Residual: " << err << std::endl;
		
		if ((condition_1) or (condition_2) or (condition_3) or (condition_4))
		{		
			is_converged = true;
			std::cout << "        b | Condition 1: " << condition_1;
			std::cout << ", Condition 2: " << condition_2;
			std::cout << ", Condition 3: " << condition_3;
			std::cout << ", Condition 4: " << condition_4;
			std::cout << std::endl;
		}
		iter++;
	}
	
	return min_err_alpha;
}

Real
ExternalRomLoopUserObject::alphaComputeGradient(Real alpha, PetscVector<Number> & soln_inp, PetscVector<Number> & del_soln)
{
	Real h = 0.00001;
	PetscVector<Number> res1(_communicator, soln_inp.size());
	PetscVector<Number> res2(_communicator, soln_inp.size());
	
	PetscVector<Number> soln1(_communicator, soln_inp.size());
	PetscVector<Number> soln2(_communicator, soln_inp.size());
	soln1.add(1, soln_inp);
	soln2.add(1, soln_inp);
	
	soln1.add(alpha, del_soln);
	soln2.add(alpha + h, del_soln);
	//updateSolutionOnBoundary(soln1, BC_map);
	//updateSolutionOnBoundary(soln2, BC_map);
	
	_fe_problem.computeResidual(soln1, *_residual, 0);
	Real err1 = _residual->l2_norm();
	
	_fe_problem.computeResidual(soln2, *_residual, 0);
	Real err2 = _residual->l2_norm();
	
	//std::cout << "          [Grad] err1: " << err1 << " | err2: " << err2 << " | h: " << h << " | grad: " << (err2 - err1) / h << std::endl;
	
	return (err2 - err1) / h;
}

Real
ExternalRomLoopUserObject::getStepSize(std::vector<Real> grad_history, Real input_step_size)
{
	if (grad_history.size() < 2)
		return input_step_size;
	else
	{
		int n = grad_history.size();
		Real current_grad = grad_history[n-1];
		Real previous_grad = grad_history[n-2];
		if (((current_grad/std::abs(current_grad)) + (previous_grad/std::abs(previous_grad))) == 0)
			return input_step_size / 2.0;
		else
			return input_step_size;
	}
}

DenseVector<Real>
ExternalRomLoopUserObject::LinearUpdateFOM(DenseMatrix<Real> * J)
{
	int max_linear_itts = _max_linear_iter;
	int linear_abs_tol =_linear_tol_ROM;
	
	bool is_converged = false;
	int l_iter = 1;
	Real err;
	Real min_res = std::numeric_limits<double>::max();
	
	DenseVector<Real> DeltaSolution(nNodes);
	_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
	DenseVector<Number> S;
	DenseVector<Number> R;
	
	DenseVector<Number> minS(S.size());
	
	std::cout << std::endl;
	while ((!is_converged) and (l_iter < max_linear_itts+1))
	{
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
		DeltaSolution.zero();
		R.scale(-1.0);
		J->lu_solve(R, DeltaSolution);
		
		Real alpha;
		if (_optimize_alpha_grad_descent)
			alpha = alphaOptimizationGradientDescent(S, DeltaSolution, _relaxation_factor);
		else if (_optimize_alpha_basic)
			alpha = getOptimalAlpha(S, DeltaSolution, _relaxation_factor);
		else
			alpha = _relaxation_factor / _alpha_opt_beta;
		S.add(alpha, DeltaSolution);
		updateSolution(S);
		
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		if (_convergence_type == 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "[Linear] Iteration " << l_iter << std::endl;
		std::cout << "        Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "        Error [" << _convergence_name << "]: " << err << std::endl;
		
		if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			minS = S;
		}
		
		l_iter++;
	}
	updateSolution(minS);
	
	return minS;
}

DenseVector<Real>
ExternalRomLoopUserObject::LinearUpdateROM(DenseMatrix<Real> * Jr)
{
	int max_linear_itts = _max_linear_iter;
	int linear_abs_tol =_linear_tol_ROM;
	
	bool is_converged = false;
	int l_iter = 1;
	Real err;
	Real min_res = std::numeric_limits<double>::max();
	
	DenseVector<Real> DeltaSolution_red(Jr->m());
	DenseVector<Real> DeltaSolution(nNodes);
	_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
	DenseVector<Number> S;
	DenseVector<Number> R;
	DenseVector<Number> Rr;
	
	DenseVector<Number> minS(S.size());
	
	S = NumericToDense(*_curr_sol);
	min_res = _residual->l2_norm();
	minS = S;
	
	std::cout << std::endl;
	while ((!is_converged) and (l_iter < max_linear_itts+1))
	{
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
		Rr = reduceVector(R);
		DeltaSolution.zero();
		R.scale(-1.0);
		Jr->lu_solve(Rr, DeltaSolution_red);
		DeltaSolution = reconstructVector(DeltaSolution_red);
		/*
		if (_scale_solution)
			if (iter == 0)
				DeltaSolution = scaleVector(DeltaSolution, _addition_scaling_vector, _multiplication_scaling_vector, false, true, true);
			else
				DeltaSolution = scaleVector(DeltaSolution, _addition_scaling_vector, _multiplication_scaling_vector, false, false, true);
		*/
		
		Real alpha;
		if (_optimize_alpha_grad_descent)
			alpha = alphaOptimizationGradientDescent(S, DeltaSolution, _relaxation_factor);
		else if (_optimize_alpha_basic)
			alpha = getOptimalAlpha(S, DeltaSolution, _relaxation_factor);
		else
			alpha = _relaxation_factor / _alpha_opt_beta;
		S.add(alpha, DeltaSolution);
		updateSolution(S);
		
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		if (_convergence_type == 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "[Linear] Iteration " << l_iter << std::endl;
		std::cout << "        Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "        Error [" << _convergence_name << "]: " << err << std::endl;
		
		if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			minS = S;
		}
		
		l_iter++;
	}
	updateSolution(minS);
	
	return minS;
}

// ******************************************************************************************************************************* //
// ****************************** LOOPS ****************************************************************************************** //
// ******************************************************************************************************************************* //

void
ExternalRomLoopUserObject::FOMLoop()
{
	// Initial print statments
	std::cout << "Full Order Model" << std::endl;
	if (_transient) {std::cout << "Transient" << std::endl; std::cout << "t_step: " << _t_step << std::endl; err_log_file << "t_step: " << _t_step << std::endl; }
	else {std::cout << "Steady-State" << std::endl;}
	// Timing
	std::string t_units = "us";
	
	// Initalize Loop Variables
	bool is_converged = false;
	iter = 0;
	Real err;
	
	// Initalize Vectors & Matrices
	DenseVector<Number> S(nNodes);  // Solution
	DenseVector<Real> R(nNodes);  // Resdiual
	SparseMatrix<Number> * J(_jac_matrix); // Jacobian
	for (int i=0; i<nNodes; i++) {S(i) = _curr_sol->el(i);}
	for (int i=0; i<nNodes; i++) {R(i) = 0;}
	DenseMatrix<Real> J_dense = SparseToDenseMatrix(J, jacobian_zeros_map);
	
	DenseVector<Number> ConvergedSolution(nNodes);  // Final converged solution
	Real min_res;
	
	while ((!is_converged) and (iter < _max_iter_rom))
	{
		// Create Dense vectors
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
		
		// Initial Log-Save
		if (_save_to_logs){
			saveVectorToLog(S, "Solution (input)", sol_log_file, true, false, iter);
			saveVectorToLog(R, "Residual (input)", res_log_file, true, false, iter);
			//saveDenseMatrixToLog(J_dense, "Jacobian (input)", jac_log_file,true, false, iter);
		}
		
		// Linear Solve
		DenseVector<Real> deltaS;
		R.scale(-1.0);
		J_dense.lu_solve(R, deltaS);
		
		
		// Update Solution		
		if (iter > 0)
		{	
			if (_use_linear_iterations)
				S = LinearUpdateFOM(&J_dense);
			else
			{
				Real alpha;
				if (_optimize_alpha_grad_descent)
					alpha = alphaOptimizationGradientDescent(S, deltaS, _relaxation_factor);
				else if (_optimize_alpha_basic)
					alpha = getOptimalAlpha(S, deltaS, _relaxation_factor);
				else
					alpha = _relaxation_factor;
				std::cout << "--- ALPHA : " << alpha << " ---" << std::endl; 
				S.add(alpha, deltaS);
			}
			
		} else {S.add(_relaxation_factor, deltaS);}
		
		
		updateSolution(S);
		
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		if (_convergence_type <= 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		else if (_convergence_type == 3)
			err = computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "Iteration " << iter << std::endl;
		std::cout << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "    Error [" << _convergence_name << "]: " << err << std::endl;
		err_log_file << "Iteration " << iter << std::endl;
		err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			err_log_file << "    Error [" << _convergence_name << "]: " << err << std::endl;
		
		//std::cout << "        (Final-" << iter << ") soln_i: [2] = " << (*_curr_sol).el(2) << " | [16] = " << (*_curr_sol).el(16)<< " | [3] = " << (*_curr_sol).el(3) << " | [13210] = " << (*_curr_sol).el(13210) << std::endl;
		
		if (iter == 0)
			min_res = _residual->l2_norm();
		else if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			ConvergedSolution = S;
		}
		
		if (err < _tolerance)
			is_converged = true;
		
		// Final Log-Save
		if (_save_to_logs){
			saveVectorToLog(deltaS, "Delta-Solution", sol_log_file, true, false, iter);
			saveVectorToLog(_curr_sol, "Solution (post-solve)", sol_log_file, true, false, iter);
			saveVectorToLog(_residual, "Residual (post-solve)", res_log_file, true, false, iter);}
		
		//for (int i=0; i<nNodes; i++)
		//	std::cout << i << ": (curr) " << (*_curr_sol).el(i) << " | (S) " << S.el(i) << " | (diff) " << (*_curr_sol).el(i) - S.el(i) << std::endl;
		
		// All Error Types
		if ((_save_all_errors) and (_convergence_type >= 0)){
			all_err_log_file << "Iteration " << iter << std::endl;
			all_err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
			all_err_log_file << "    RMSE: " << computeRMSE(_expected_soln, S) << std::endl;
			all_err_log_file << "    MaxError: " << computeMaxError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;
			all_err_log_file << "    MaxPercentError: " << computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;}
		
		// Save solution vector
		if (_save_final_solution)
			saveVectorToLog(_curr_sol, "Solution", final_sol_log_file, true, false, iter);
				
		iter++;
		
		if (!is_converged)
		{
			_nl_sys.computeJacobian(*J);
			J_dense = SparseToDenseMatrix(J, jacobian_zeros_map);
		}
				
		
		if (_save_timing)
		{
			
		}
	}
	
	if (!is_converged)
	{
		_t = _t - _dt;
		_dt = _dt / 2;
		std::cout << "ERROR: DID NOT CONVERGE" << std::endl;
	} 
	
	saveVectorToLog(ConvergedSolution, "Solution", conv_sol_log_file, true, false, _t_step);
	updateSolution(ConvergedSolution);
	std::cout << "Minimum Residual: " << min_res << std::endl;
	
}

void
ExternalRomLoopUserObject::ROMLoop()
{
	// Initial print statments
	std::cout << "Reduced Order Model" << std::endl;
	if (_transient) {std::cout << "Transient" << std::endl; std::cout << "t_step: " << _t_step << std::endl; err_log_file << "t_step: " << _t_step << std::endl; }
	else {std::cout << "Steady-State" << std::endl;}
	// Timing
	std::string t_units = "us";
	
	// Initalize Loop Variables
	bool is_converged = false;
	iter = 0;
	Real err;
	
	// Initalize Vectors & Matrices
	DenseVector<Number> S(nNodes);  // Solution
	DenseVector<Number> R(nNodes);  // Resdiual
	SparseMatrix<Number> * J(_jac_matrix); // Jacobian
	
	DenseVector<Number> Sr;
	DenseVector<Number> Rr;
	DenseMatrix<Number> Jr;
	
	DenseVector<Number> ConvergedSolution(nNodes);  // Final converged solution
	Real min_res;
	
	//if (_scale_solution)
	//	S = scaleVector(S, _addition_scaling_vector, _multiplication_scaling_vector, false, true, false);
	
	while ((!is_converged) and (iter < _max_iter_rom))
	{
		int loop_start = currentTime(t_units);
		// Create Dense vectors
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
		
		// Initial Reduction
		int reduction_start = currentTime(t_units);
		//Sr = reduceVector(S);
		Rr = reduceVector(R);
		Jr = reduceMatrix(J, jacobian_zeros_map);
		int reduction_end = currentTime(t_units);
		
		// Initial Log-Save
		if (_save_to_logs){
			saveVectorToLog(S, "Solution (input)", sol_log_file, true, false, iter);
			//saveVectorToLog(Sr, "Reduced Solution (input)", sol_log_file, true, false, iter);
			saveVectorToLog(R, "Residual (input)", res_log_file, true, false, iter);
			saveVectorToLog(Rr, "Reduced Residual (input)", res_log_file, true, false, iter);
			//saveSparseMatrixToLog(J, "Jacobian (input)", jac_log_file, true, false, iter);
			//saveDenseMatrixToLog(Jr, "Jacobian (input)", jac_log_file,true, false, iter);
		}
		
		// Linear Solve
		int lu_solve_start = currentTime(t_units);
		DenseVector<Number> deltaSr;
		Rr.scale(-1.0);
		Jr.lu_solve(Rr, deltaSr);
		int lu_solve_end = currentTime(t_units);
		
		// Update Solution
		int soln_update_start = currentTime(t_units);
		DenseVector<Number> deltaS = reconstructVector(deltaSr);
		if (_scale_solution)
			if (iter == 0)
				deltaS = scaleVector(deltaS, _addition_scaling_vector, _multiplication_scaling_vector, false, true, true);
			else
				deltaS = scaleVector(deltaS, _addition_scaling_vector, _multiplication_scaling_vector, false, false, true);
		
		/*
		Real alpha;
		if (iter > 0)
		{
			if (_optimize_alpha_grad_descent)
				alpha = alphaOptimizationGradientDescent(S, deltaS, _relaxation_factor);
			else if (_optimize_alpha_basic)
				alpha = getOptimalAlpha(S, deltaS, _relaxation_factor);
			else
				alpha = _relaxation_factor;
		} else {alpha = _relaxation_factor;}
		*/
		
		// Update Solution		
		if (iter > 0)
		{	
			if (_use_linear_iterations)
				S = LinearUpdateROM(&Jr);
			else
			{
				Real alpha;
				if (_optimize_alpha_grad_descent)
					alpha = alphaOptimizationGradientDescent(S, deltaS, _relaxation_factor);
				else if (_optimize_alpha_basic)
					alpha = getOptimalAlpha(S, deltaS, _relaxation_factor);
				else
					alpha = _relaxation_factor;
				std::cout << "--- ALPHA : " << alpha << " ---" << std::endl; 
				S.add(alpha, deltaS);
			}
			
		} else {S.add(_relaxation_factor, deltaS);}
		
		updateSolution(S);
		int soln_update_end = currentTime(t_units);
		
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		if (_convergence_type <= 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		else if (_convergence_type == 3)
			err = computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "Iteration " << iter << std::endl;
		std::cout << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "    Error [" << _convergence_name << "]: " << err << std::endl;
		err_log_file << "Iteration " << iter << std::endl;
		err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			err_log_file << "    Error [" << _convergence_name << "]: " << err << std::endl;
		
		if (iter == 0)
			min_res = _residual->l2_norm();
		else if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			ConvergedSolution = S;
		}
		
		if (err < _tolerance)
			is_converged = true;
		
		// Final Log-Save
		if (_save_to_logs){
			saveVectorToLog(deltaSr, "Reduced Delta-Solution", sol_log_file, true, false, iter);
			saveVectorToLog(reconstructVector(deltaSr), "Delta-Solution", sol_log_file, true, false, iter);
			if (_scale_solution)
				saveVectorToLog(deltaS, "Delta-Solution (Scaled)", sol_log_file, true, false, iter);
			saveVectorToLog(_curr_sol, "Solution (post-solve)", sol_log_file, true, false, iter);
			saveVectorToLog(_residual, "Residual (post-solve)", res_log_file, true, false, iter);}
		
		// All Error Types
		if ((_save_all_errors) and (_convergence_type >= 0)){
			all_err_log_file << "Iteration " << iter << std::endl;
			all_err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
			all_err_log_file << "    RMSE: " << computeRMSE(_expected_soln, S) << std::endl;
			all_err_log_file << "    MaxError: " << computeMaxError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;
			all_err_log_file << "    MaxPercentError: " << computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;}
		
		// Save solution vector
		if (_save_final_solution)
			saveVectorToLog(S, "Solution", final_sol_log_file, true, false, iter);
		
		iter++;
		
		int jac_update_start = currentTime(t_units);
		if (!is_converged)
		{
			_nl_sys.computeJacobian(*J);
		}
		int jac_update_end = currentTime(t_units);
		
		int loop_end = currentTime(t_units);
		
		if (_save_timing)
		{
			int loop_time        = loop_end - loop_start;
			int reduction_time   = reduction_end - reduction_start;
			int lu_solve_time    = lu_solve_end - lu_solve_start;
			int soln_update_time = soln_update_end - soln_update_start;
			int jac_update_time  = jac_update_end - jac_update_start;
			
			saveStringToTimeLog("Reduction", reduction_time, t_units);
			saveStringToTimeLog("LU_Solve", lu_solve_time, t_units);
			saveStringToTimeLog("Solution_Update", soln_update_time, t_units);
			saveStringToTimeLog("Jacobian_Update", jac_update_time, t_units);
		}
	}
	
	saveVectorToLog(ConvergedSolution, "Solution", conv_sol_log_file, true, false, _t_step);
	
	if (!is_converged)
		std::cout << "ERROR: DID NOT CONVERGE" << std::endl;
}

void
ExternalRomLoopUserObject::FOMLoopDEBUG()
{
	// Initial print statments
	std::cout << "\nFull Order Model (DEBUG)" << std::endl;
	if (_transient) {std::cout << "Transient" << std::endl; std::cout << "t_step: " << _t_step << std::endl; err_log_file << "t_step: " << _t_step << std::endl; }
	else {std::cout << "Steady-State" << std::endl;}
	// Timing
	std::string t_units = "us";
	
	// Initalize Loop Variables
	bool is_converged = false;
	iter = 0;
	Real err;
	
	// Initalize Vectors & Matrices
	DenseVector<Number> S(nNodes);  // Solution
	DenseVector<Real> R(nNodes);  // Resdiual
	SparseMatrix<Number> * J(_jac_matrix); // Jacobian
	for (int i=0; i<nNodes; i++) {S(i) = _curr_sol->el(i);}
	for (int i=0; i<nNodes; i++) {R(i) = 0;}
	std::cout << "\nHERE1" << std::endl;
	DenseMatrix<Real> J_dense = SparseToDenseMatrix(J, jacobian_zeros_map);
	std::cout << "HERE2" << std::endl;
	
	DenseVector<Number> ConvergedSolution(nNodes);  // Final converged solution
	Real min_res;
	
	
	while ((!is_converged) and (iter < _max_iter_rom))
	{
		// Create Dense vectors
		//S = NumericToDense(*_curr_sol);
		//R = NumericToDense(*_residual);
		
		// Initial Log-Save
		/*
		if (_save_to_logs){
			saveVectorToLog(S, "Solution (input)", sol_log_file, true, false, iter);
			saveVectorToLog(R, "Residual (input)", res_log_file, true, false, iter);
			//saveDenseMatrixToLog(J_dense, "Jacobian (input)", jac_log_file,true, false, iter);
		}
		*/
		
		// Linear Solve
		DenseVector<Real> deltaS;
		R.scale(-1.0);
		std::cout << "\nHere1" << std::endl;
		(&J_dense)->lu_solve(R, deltaS);
		std::cout << "\nHere2" << std::endl;
		//Real tol2 = std::pow(10, -10);
		//int m_its = 10;
		//_residual->scale(-1.0);
		
		/*
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		std::cout << "initial: "  << " | R = " << _residual->l2_norm() << std::endl;
		linear_solver.clear();
		for (int i=0; i<m_its; i++)
		{
			linear_solver.solve(*J, *_curr_sol, *_residual, tol2, _max_linear_iter);
			S = NumericToDense(*_curr_sol);
			updateSolution(S);
			_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
			std::cout << "i: " << i << " | R = " << _residual->l2_norm() << std::endl;
		}
		*/
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
		std::cout << "\nHere3" << std::endl;
		// Update Solution		
		/*
		if (iter > 0)
		{	
			if (_use_linear_iterations)
				S = LinearUpdateFOM(&J_dense);
			else
			{
				Real alpha;
				if (_optimize_alpha_grad_descent)
					alpha = alphaOptimizationGradientDescent(S, deltaS, _relaxation_factor);
				else if (_optimize_alpha_basic)
					alpha = getOptimalAlpha(S, deltaS, _relaxation_factor);
				else
					alpha = _relaxation_factor;
				std::cout << "--- ALPHA : " << alpha << " ---" << std::endl; 
				S.add(alpha, deltaS);
			}
			
		} else {S.add(_relaxation_factor, deltaS);}
		*/
		updateSolution(S);
		std::cout << "\nHere4" << std::endl;
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		if (_convergence_type <= 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		else if (_convergence_type == 3)
			err = computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "Iteration " << iter << std::endl;
		std::cout << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "    Error [" << _convergence_name << "]: " << err << std::endl;
		err_log_file << "Iteration " << iter << std::endl;
		err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			err_log_file << "    Error [" << _convergence_name << "]: " << err << std::endl;
		
		//std::cout << "        (Final-" << iter << ") soln_i: [2] = " << (*_curr_sol).el(2) << " | [16] = " << (*_curr_sol).el(16)<< " | [3] = " << (*_curr_sol).el(3) << " | [13210] = " << (*_curr_sol).el(13210) << std::endl;
		
		if (iter == 0)
			min_res = _residual->l2_norm();
		else if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			ConvergedSolution = S;
		}
		
		if (err < _tolerance)
			is_converged = true;
		
		// Final Log-Save
		if (_save_to_logs){
			//saveVectorToLog(deltaS, "Delta-Solution", sol_log_file, true, false, iter);
			saveVectorToLog(_curr_sol, "Solution (post-solve)", sol_log_file, true, false, iter);
			saveVectorToLog(_residual, "Residual (post-solve)", res_log_file, true, false, iter);}
		
		//for (int i=0; i<nNodes; i++)
		//	std::cout << i << ": (curr) " << (*_curr_sol).el(i) << " | (S) " << S.el(i) << " | (diff) " << (*_curr_sol).el(i) - S.el(i) << std::endl;
		
		// All Error Types
		if ((_save_all_errors) and (_convergence_type >= 0)){
			all_err_log_file << "Iteration " << iter << std::endl;
			all_err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
			all_err_log_file << "    RMSE: " << computeRMSE(_expected_soln, S) << std::endl;
			all_err_log_file << "    MaxError: " << computeMaxError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;
			all_err_log_file << "    MaxPercentError: " << computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;}
		
		// Save solution vector
		if (_save_final_solution)
			saveVectorToLog(_curr_sol, "Solution", final_sol_log_file, true, false, iter);
				
		iter++;
		
		if (!is_converged)
		{
			std::cout << "\nHERE5" << std::endl;
			_nl_sys.computeJacobian(*J);
			std::cout << "HERE6" << std::endl;
			J_dense = SparseToDenseMatrix(J, jacobian_zeros_map);
			std::cout << "HERE7" << std::endl;
		}
				
		
		if (_save_timing)
		{
			
		}
	}
	
	if (!is_converged)
		std::cout << "ERROR: DID NOT CONVERGE" << std::endl;
	
	saveVectorToLog(ConvergedSolution, "Solution", conv_sol_log_file, true, false, _t_step);
	updateSolution(ConvergedSolution);
	std::cout << "Minimum Residual: " << min_res << std::endl;
	
}

void
ExternalRomLoopUserObject::ROMLoopDEBUG()
{
	// Initial print statments
	std::cout << "Reduced Order Model (DEBUG)" << std::endl;
	if (_transient) {std::cout << "Transient" << std::endl; std::cout << "t_step: " << _t_step << std::endl; err_log_file << "t_step: " << _t_step << std::endl; }
	else {std::cout << "Steady-State" << std::endl;}
	// Timing
	std::string t_units = "us";
	
	// Initalize Loop Variables
	bool is_converged = false;
	iter = 0;
	Real err;
	
	// Initalize Vectors & Matrices
	DenseVector<Number> S(nNodes);  // Solution
	DenseVector<Number> R(nNodes);  // Resdiual
	SparseMatrix<Number> * J(_jac_matrix); // Jacobian
		
	DenseVector<Number> Sr;
	DenseVector<Number> Rr;
	DenseMatrix<Number> Jr;
	
	DenseVector<Number> ConvergedSolution(nNodes);  // Final converged solution
	Real min_res;
	
	if (_scale_solution)
	{
		S.zero();
		S = scaleVector(S, _addition_scaling_vector, _multiplication_scaling_vector, false, true, false);
		updateSolution(S);
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
	}
		
	//if (_scale_solution)
	//	S = scaleVector(S, _addition_scaling_vector, _multiplication_scaling_vector, false, true, false);
	
	while ((!is_converged) and (iter < _max_iter_rom))
	{
		int loop_start = currentTime(t_units);
		// Create Dense vectors
		S = NumericToDense(*_curr_sol);
		R = NumericToDense(*_residual);
				
		// Initial Reduction
		int reduction_start = currentTime(t_units);
		Sr = reduceVector(S);
		Rr = reduceVector(R);
		int jac_reduction_start = currentTime(t_units);
		Jr = reduceMatrix(J, jacobian_zeros_map);
		int jac_reduction_end = currentTime(t_units);
		int reduction_end = currentTime(t_units);
				
		// Initial Log-Save
		if (_save_to_logs){
			saveVectorToLog(S, "Solution (input)", sol_log_file, true, false, iter);
			saveVectorToLog(Sr, "Reduced Solution (input)", sol_log_file, true, false, iter);
			saveVectorToLog(R, "Residual (input)", res_log_file, true, false, iter);
			saveVectorToLog(Rr, "Reduced Residual (input)", res_log_file, true, false, iter);
			saveSparseMatrixToLog(J, "Jacobian (input)", jac_log_file, true, false, iter);
			saveDenseMatrixToLog(Jr, "Jacobian (input)", jac_log_file,true, false, iter);}
		
		// Linear Solve
		int lu_solve_start = currentTime(t_units);
		DenseVector<Number> deltaSr;
		Rr.scale(-1.0);
		Jr.lu_solve(Rr, deltaSr);
		int lu_solve_end = currentTime(t_units);
		
		// Update Solution
		int soln_update_start = currentTime(t_units);
		DenseVector<Number> deltaS = reconstructVector(deltaSr);
		/*
		if (_scale_solution)
			if (iter == 0)
				deltaS = scaleVector(deltaS, _addition_scaling_vector, _multiplication_scaling_vector, false, true, true);
			else
				deltaS = scaleVector(deltaS, _addition_scaling_vector, _multiplication_scaling_vector, false, false, true);
		
		debug_log_file << "SOLUTION || Iter : " << iter << std::endl;
		debug_log_file << "dof : S(i) | delS(i) | (S+delS)(i) | Strue(i)" << std::endl;
		for (int i=2379; i<3013; i++)
		{
			debug_log_file << i << " : " << S.el(i) << " | " << deltaS.el(i) << " | " << S.el(i) + deltaS.el(i) << " | " << _expected_soln.el(i) << std::endl;
		}
		*/
		
		// Update Solution		
		if (iter > 0)
		{	
			if (_use_linear_iterations)
				S = LinearUpdateROM(&Jr);
			else
			{
				Real alpha;
				if (_optimize_alpha_grad_descent)
					alpha = alphaOptimizationGradientDescent(S, deltaS, _relaxation_factor);
				else if (_optimize_alpha_basic)
					alpha = getOptimalAlpha(S, deltaS, _relaxation_factor);
				else
					alpha = _relaxation_factor;
				std::cout << "--- ALPHA : " << alpha << " ---" << std::endl; 
				S.add(alpha, deltaS);
			}
			
		} else {S.add(_relaxation_factor, deltaS);}
		
		updateSolution(S);
		int soln_update_end = currentTime(t_units);
		
		// Compute Errors
		_fe_problem.computeResidual(*_curr_sol, *_residual, 0);
		
		debug_log_file << "RESDIUAL || Iter : " << iter << std::endl;
		debug_log_file << "dof : R_t(i) | R_t+1(i)" << std::endl;
		for (int i=2379; i<3013; i++)
		{
			debug_log_file << i << " : " << R.el(i) << " | " << _residual->el(i) << std::endl;
		}
		
		if (_convergence_type <= 0)
			err = _residual->l2_norm();
		else if (_convergence_type == 1)
			err = computeRMSE(_expected_soln, S);
		else if (_convergence_type == 2)
			err = computeMaxError(_expected_soln, NumericToDense(*_curr_sol));
		else if (_convergence_type == 3)
			err = computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol));
		
		std::cout << "Iteration " << iter << std::endl;
		std::cout << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			std::cout << "    Error [" << _convergence_name << "]: " << err << std::endl;
		err_log_file << "Iteration " << iter << std::endl;
		err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
		if (_convergence_type > 0)
			err_log_file << "    Error [" << _convergence_name << "]: " << err << std::endl;
		
		if (iter == 0)
			min_res = _residual->l2_norm();
		else if (_residual->l2_norm() < min_res)
		{
			min_res = _residual->l2_norm();
			ConvergedSolution = S;
		}
		
		if (err < _tolerance)
			is_converged = true;
		
		// Final Log-Save
		if (_save_to_logs){
			saveVectorToLog(deltaSr, "Reduced Delta-Solution", sol_log_file, true, false, iter);
			saveVectorToLog(reconstructVector(deltaSr), "Delta-Solution", sol_log_file, true, false, iter);
			if (_scale_solution)
				saveVectorToLog(deltaS, "Delta-Solution (Scaled)", sol_log_file, true, false, iter);
			saveVectorToLog(_curr_sol, "Solution (post-solve)", sol_log_file, true, false, iter);
			saveVectorToLog(_residual, "Residual (post-solve)", res_log_file, true, false, iter);}
		
		// All Error Types
		if ((_save_all_errors) and (_convergence_type >= 0)){
			all_err_log_file << "Iteration " << iter << std::endl;
			all_err_log_file << "    Residual: " << _residual->l2_norm() << std::endl;
			all_err_log_file << "    RMSE: " << computeRMSE(_expected_soln, S) << std::endl;
			all_err_log_file << "    MaxError: " << computeMaxError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;
			all_err_log_file << "    MaxPercentError: " << computeMaxPercentError(_expected_soln, NumericToDense(*_curr_sol)) << std::endl;}
		
		// Save solution vector
		if (_save_final_solution)
			saveVectorToLog(S, "Solution", final_sol_log_file, true, false, iter);
		
		iter++;
		
		int jac_update_start = currentTime(t_units);
		if (!is_converged)
		{
			_nl_sys.computeJacobian(*J);
		}
		int jac_update_end = currentTime(t_units);
		
		int loop_end = currentTime(t_units);
		
		if (_save_timing)
		{
			int loop_time        = loop_end - loop_start;
			int reduction_time   = reduction_end - reduction_start;
			int jac_reduction_time   = jac_reduction_end - jac_reduction_start;
			int lu_solve_time    = lu_solve_end - lu_solve_start;
			int soln_update_time = soln_update_end - soln_update_start;
			int jac_update_time  = jac_update_end - jac_update_start;
			
			saveStringToTimeLog("Reduction", reduction_time, t_units);
			saveStringToTimeLog("Jacobian_Reduction", jac_reduction_time, t_units);
			saveStringToTimeLog("LU_Solve", lu_solve_time, t_units);
			saveStringToTimeLog("Solution_Update", soln_update_time, t_units);
			saveStringToTimeLog("Jacobian_Update", jac_update_time, t_units);
		}
	}
	
	saveVectorToLog(ConvergedSolution, "Solution", conv_sol_log_file, true, false, _t_step);
	
	if (!is_converged)
		std::cout << "ERROR: DID NOT CONVERGE" << std::endl;
}


// ******************************************************************************************************************************* //