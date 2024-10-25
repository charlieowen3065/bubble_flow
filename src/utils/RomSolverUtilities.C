#include "RomSolverUtilities.h"
#include "DelimitedFileReader.h"
#include "petscmat.h" 

RomSolverUtilities::RomSolverUtilities()
{
	
}


// ************************* Getters & Setters ************************* //
void 
RomSolverUtilities::setLogFiles(std::vector<std::string> filenames, bool save, bool timer)
{
	gen_log_file.open(filenames[0]); 
	sol_log_file.open(filenames[1]);
	res_log_file.open(filenames[2]);
	jac_log_file.open(filenames[3]);
	
	time_log_file.open(filenames[4]);
	
	_save_to_logs = save;
	_save_timing = timer; 
}

// ************************* Setup Functions ************************* //
void
RomSolverUtilities::loadQMatrix(std::string filepath, int nRanks_input, const libMesh::Parallel::Communicator * comm)
{	
	std::cout << "loadQMatrix()" << std::endl;
	
	MooseUtils::DelimitedFileReader Q_csv_reader(filepath, comm);
	Q_csv_reader.read(); const std::vector<std::vector<double>> & Q_inp = Q_csv_reader.getData();
	nNodes = Q_inp[0].size();
	if (nRanks_input == -1) {nRanks = Q_inp.size();}
	else {nRanks = nRanks_input;}
	
	std::cout << std::endl << "nNodes: " << nNodes << std::endl << "nRanks: " << nRanks << std::endl;

	Q_dense.resize(nNodes, nRanks);
	for (int i = 0; i < nNodes; i++)
		for (int j = 0; j < nRanks; j++)
			Q_dense(i, j) = Q_inp[j][i];
}

DenseVector<Real>
RomSolverUtilities::loadVector(std::string filepath, const libMesh::Parallel::Communicator * comm)
{
	MooseUtils::DelimitedFileReader V_csv_reader(filepath, comm);
	V_csv_reader.read(); const std::vector<std::vector<double>> & V_inp = V_csv_reader.getData();
	
	int N = V_inp[0].size(); DenseVector<Real> V(N);
	for (int i=0; i<N; i++)
		V(i) = V_inp[0][i];
	
	if (_save_to_logs)
		saveVectorToLog(V, filepath, gen_log_file, true, false, 0);
	
	return V;	
}

// ************************* Log-File Functions ************************* //
void
RomSolverUtilities::saveStdMatrixToLog(std::vector<std::vector<Real>> matrix, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter)
{
	int n;;
	if (!transpose) {n = matrix[0].size();} 
	else {n = matrix.size();}
	
	if (include_iter) {log_file << header << " (iter = " << iter << ")" << std::endl;}
	else {log_file << header << std::endl;}
	
	log_file << "Row | Column Values ..." << std::endl;

	std::string lines;	
	for (int i=0; i<n; i++)
	{
		lines += std::to_string(i) + ": ";
		for (auto val : matrix[i])
			lines += std::to_string(val) + " ";
		lines += "\n";
	}
	log_file << lines << std::endl;
}

void
RomSolverUtilities::saveSparseMatrixToLog(SparseMatrix<Real> * matrix, std::string header, std::ofstream & log_file, bool include_iter, bool sparse, int iter)
{
	if (include_iter) {log_file << header << " (iter = " << iter << ")" << std::endl;}
	else {log_file << header << std::endl;}
	matrix->print(log_file, sparse);
	log_file << std::endl;
}

void
RomSolverUtilities::saveDenseMatrixToLog(DenseMatrix<Real> & matrix, std::string header, std::ofstream & log_file, bool include_iter, bool sparse, int iter)
{
	if (include_iter) {log_file << header << " (iter = " << iter << ")" << std::endl;}
	else {log_file << header << std::endl;}
	
	log_file << "Row | Column Values ..." << std::endl;
	
	int n = matrix.n(/*Columns*/); int m = matrix.m(/*Rows*/);  
	
	std::string lines;
	for (int i=0; i<m; i++)
	{
		lines += std::to_string(i) + ": ";
		for (int j=0; j<n; j++)
		{
			Real val = matrix(i,j);
			lines += std::to_string(val) + " ";
		}
		lines += "\n";
	}
	log_file << lines << std::endl;
}

void
RomSolverUtilities::saveVectorToLog(DenseVector<Real> vector, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter)
{	
	if (include_iter) {log_file << header << " (iter = " << iter << ")" << std::endl;}
	else {log_file << header << std::endl;}
	
	log_file << "Row | Column Values ..." << std::endl;
	
	std::string lines;
	int nrows = vector.size();
	for (int i=0; i<nrows; i++)
	{
		lines += std::to_string(i) + ": ";
		lines += std::to_string(vector.el(i));
		lines += "\n";
	}
	log_file << lines << std::endl;
}

void
RomSolverUtilities::saveVectorToLog(const NumericVector<Real> * vector, std::string header, std::ofstream & log_file, bool include_iter, bool transpose, int iter)
{
	if (include_iter) {log_file << header << " (iter = " << iter << ")" << std::endl;}
	else {log_file << header << std::endl;}
	
	log_file << "Row | Column Values ..." << std::endl;
	int nrows = vector->size();
	
	std::string lines;
	for (int i=0; i<nrows; i++)
	{
		lines += std::to_string(i) + ": ";
		lines += std::to_string(vector->el(i));
		lines += "\n";
	}
	log_file << lines << std::endl;
}

// ************************* Timing Functions ************************* //
int
RomSolverUtilities::currentTime(std::string units)
{
	int now;
	if ((units == "s") or (units == "sec") or (units == "seconds"))
		now = currentTime_s();
	else if ((units == "ms") or (units == "millisec") or (units == "milliseconds"))
		now = currentTime_ms();
	else if ((units == "us") or (units == "microsec") or (units == "microseconds"))
		now = currentTime_us();
	else
		mooseError("Please use proper units for time [s, ms, us]");

	return now;
}

// ************************* Error Functions ************************* //
Real
RomSolverUtilities::computeRMSE(DenseVector<Real> vec1, DenseVector<Real> vec2)
{
	int N; int N1 = vec1.size(); int N2 = vec2.size();
	if (N1 == N2) {N = N1;}
	else {
		std::cout << "N1: " << N1 << std::endl;
		std::cout << "N2: " << N2 << std::endl;
		//for (int i=0; ((i<N1) and (i<N2)); i++)
		//	std::cout << i << ": (1) " << vec1.el(i) << " | (2) " << vec2.el(i) << std::endl; 
		for (int i=0; i<N1; i++)
			std::cout << i << ": (1) " << vec1.el(i) << std::endl;
		for (int i=0; i<N2; i++)
			std::cout << i << ": (2) " << vec2.el(i) << std::endl;
		mooseError("Error: Vector sizes not compatable (computeRMSE)");
		}
	
	Real sum = 0.0;
	for (int i=0; i<N; i++)
		sum += std::pow(vec1.el(i)-vec2.el(i), 2);
	Real rmse = std::sqrt(sum / (Real)N);
	return rmse;	
}

Real
RomSolverUtilities::computeMaxError(DenseVector<Real> vec1, DenseVector<Real> vec2)
{
	int N; int N1 = vec1.size(); int N2 = vec2.size();
	if (N1 == N2) {N = N1;}
	else {mooseError("Error: Vector sizes not compatable (computeMaxError)");}
	
	Real max_error = 0.0;
	for (int i=0; i<N; i++)
		if ((std::abs(vec1.el(i) - vec2.el(i))) > max_error)
			max_error = std::abs(vec1.el(i) - vec2.el(i));
	return max_error;	
}

Real
RomSolverUtilities::computeMaxPercentError(DenseVector<Real> actual, DenseVector<Real> theoretical)
{
	int N; int N1 = actual.size(); int N2 = theoretical.size();
	
	if (N1 == N2) {N = N1;}
	else {mooseError("Error: Vector sizes not compatable (computeMaxError)");}
	
	std::cout << std::endl << std::endl;
	
	Real tol = 1*std::pow(10, -6);
	
	Real p_err;
	Real max_percent_error = 0.0;
	for (int i=0; i<N; i++)
	{
		Real t = theoretical.el(i);
		Real a = actual.el(i);
		
		if ((std::abs(t)<tol) or (std::abs(a)<tol))
		{
			p_err = std::abs(t-a) * 100;
			/*
			if ((std::abs(t) > 0.001) or (std::abs(a) > 0.001))
				p_err = std::abs(t-a) * 100;
			else
				p_err = 0;
			*/
		}
		else
			p_err = 100.0*(std::abs(a - t)/t);
		
		if (p_err > max_percent_error)
		{
			max_percent_error = p_err;
		}
	}
	return max_percent_error;	
}


// ************************* Mapping Functions ************************* //
std::vector<std::vector<Real>>
RomSolverUtilities::getJacobianRowVectorNonzeroIndices(SparseMatrix<Real> * jacobian)
{
	int n = jacobian->n(/*Columns*/); int m = jacobian->m(/*Rows*/);
	std::vector<std::vector<Real>> row_vector_indicies;
	for (numeric_index_type i=0; i < m; i++)
	{
		std::vector<Real> row_i_indicies;
		for (numeric_index_type j=0; j < n; j++)
		{
			if ((*jacobian)(i, j) != 0)
				row_i_indicies.push_back(j);
		}
		row_vector_indicies.push_back(row_i_indicies);
	}
	return row_vector_indicies;
}

std::map<int, std::vector<int>>
RomSolverUtilities::getJacobianZerosMap(SparseMatrix<Real> * jacobian)
{
	int n = jacobian->n(/*Columns*/); int m = jacobian->m(/*Rows*/);
	std::map<int, std::vector<int>> zeros_map;
	for (int i=0; i < m; i++)
	{
		std::vector<int> row_i_indicies;
		for (int j=0; j < n; j++)
		{
			if ((*jacobian)(i, j) != 0)
				row_i_indicies.push_back(j);
		}
		zeros_map[i] = row_i_indicies;
	}
	return zeros_map;
}

std::map<int, std::pair<std::string, Real>>
RomSolverUtilities::getBCMap(SparseMatrix<Number> * jacobian, std::map<std::string, std::vector<dof_id_type>> dof_map, DenseVector<Real> initial_solution)
{
	// Get Variables Names from DOF Map
	std::vector<std::string> var_names;
	for(std::map<std::string, std::vector<dof_id_type>>::iterator it = dof_map.begin(); it != dof_map.end(); ++it)
		var_names.push_back(it->first);
	// Get BC nodes and corresponding BC value
	std::map<int, std::pair<std::string, Real>> BC_map;
	for (auto v : var_names)
	{
		std::vector<dof_id_type> dof_ids = dof_map[v];
		for (auto dof_i : dof_ids)
		{
			Real j_ii = (*jacobian)(dof_i,dof_i);
			if (j_ii == 1)
			{
				std::pair<std::string, Real> pair_i = std::make_pair(v, initial_solution.el(dof_i));
				BC_map[dof_i] = pair_i;
			}
		}
	}
	
	return BC_map;
}

std::map<std::string, std::vector<dof_id_type>>
RomSolverUtilities::getDOFMapping_Old(std::vector<MooseVariableFieldBase*> & variables, MooseMesh & mesh)
{
	bool is_scalar;
	std::vector<std::string> comp_suffixs{ "_x", "_y", "_z" };
	std::vector<dof_id_type> temp_vec;
	
	for (auto v : variables)
	{
		std::string v_name = v->name(); 
		std::cout << "VARIABLE: " << v_name << std::endl;
		
		int num_componets = mesh.nodePtr(0)->n_comp(0, v->number());
		int num_nodes = mesh.nNodes();
		
		if (num_componets == 1) {is_scalar = true;}
		else {is_scalar = false;}
		
		if (is_scalar)
		{
			temp_vec.clear();
			for (auto node_idx : make_range(num_nodes))
			{
				const Node * node_i = mesh.nodePtr(node_idx);
				dof_id_type dof_i = node_i->dof_number(0, v->number(), 0);
				temp_vec.push_back(dof_i);
			}
			dof_map_old[v_name] = temp_vec;
		}
		else
		{
			for (auto comp_number : make_range(num_componets))
			{
				temp_vec.clear();
				for (auto node_idx : make_range(num_nodes))
				{
					const Node * node_i = mesh.nodePtr(node_idx);
					dof_id_type dof_i = node_i->dof_number(0, v->number(), comp_number);
					temp_vec.push_back(dof_i);
				}
				dof_map_old[v_name + comp_suffixs[comp_number]] = temp_vec;
			}
		}
		
	}
	return dof_map_old;
}

std::map<dof_id_type, std::string>
RomSolverUtilities::getDOFMapping(std::vector<MooseVariableFieldBase*> & variables, MooseMesh & mesh)
{
	std::cout << "********** INSIDE DOF MAPPING **********" << std::endl << std::endl;
	
	bool is_scalar;
	std::vector<std::string> comp_suffixs{ "_x", "_y", "_z" };
	for (auto v : variables)
	{
		std::string v_name = v->name(); 
		std::cout << "VARIABLE: " << v_name << std::endl;
		
		int num_componets = mesh.nodePtr(0)->n_comp(0, v->number());
		int num_nodes = mesh.nNodes();
		
		if (num_componets == 1) {is_scalar = true;}
		else {is_scalar = false;}
		
		for (int n=0; n<num_componets; n++)
		{
			std::string v_name_inp;
			if (!is_scalar)
				v_name_inp = v_name + comp_suffixs[n];
			else
				v_name_inp = v_name;
			for (auto node_idx : make_range(num_nodes))
			{
				const Node * node_i = mesh.nodePtr(node_idx);
				dof_id_type dof_i = node_i->dof_number(0, v->number(), n);
				dof_map[dof_i] = v_name_inp;
			}
		}
	}
	
	//for(std::map<dof_id_type, std::string>::iterator it = dof_map.begin(); it != dof_map.end(); ++it)
	//	std::cout << "DOF: " << it->first << ", Variable: " << it->second << std::endl;
		
	std::cout << std::endl << "********** END OF DOF MAPPING **********" << std::endl;
	return dof_map;
	
}

// ************************* Type-Transfer Functions ************************* //
DenseMatrix<Number>
RomSolverUtilities::SparseToDenseMatrix(SparseMatrix<Number> * S)
{
	int n = S->n(/*Columns*/); int m = S->m(/*Rows*/);  
	std::vector<std::vector<Real>> D_mat;
	for (int i=0; i<m; i++)
	{
		std::vector<Real> row_i;
		for (int j=0; j<n; j++)
			row_i.push_back((*S)(i,j));
		D_mat.push_back(row_i);
	}
	DenseMatrix<Number> D(m, n);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			D(i, j) = D_mat[i][j];
	return D;
}

DenseMatrix<Number> 
RomSolverUtilities::SparseToDenseMatrix(SparseMatrix<Number> * S, std::map<int, std::vector<int>> matrix_zeros_map)
{
	int n = S->n(/*Columns*/); int m = S->m(/*Rows*/); 
	DenseMatrix<Number> D(m, n);
	for (int i=0; i<m; i++)
		for (auto p_i : matrix_zeros_map[i])
			D(i, p_i) = (*S)(i, p_i);
	return D;
}

DenseVector<Number>
RomSolverUtilities::NumericToDense(NumericVector<Real> & N)
{
	int s = N.size();
	DenseVector<Real> D(s);
	for (int i=0; i<s; i++)
		D(i) = N.el(i);
	return D;
}


// ************************* Reduction Functions ************************* //
DenseMatrix<Number> 
RomSolverUtilities::computeJQ(SparseMatrix<Real> * J, std::map<int, std::vector<int>> row_vector_indicies)
{	
	DenseMatrix<Real> JQ(nNodes, nRanks);
	for (int i=0; i<nNodes; i++)
	{
		for (int j=0; j<nRanks; j++)
		{
			float sum = 0;
			for (auto p_i : row_vector_indicies[i])
				sum += (*J)(i, p_i) * Q_dense.el(p_i, j);
			JQ(i, j) = sum;
		}
	}
	return JQ;
}

DenseMatrix<Number>
RomSolverUtilities::reduceMatrix(SparseMatrix<Real> * M, std::map<int, std::vector<int>> row_vector_indicies)
{
	DenseMatrix<Real> M_dense(nNodes, nRanks);
	M_dense = computeJQ(M, row_vector_indicies);
	M_dense.left_multiply_transpose(Q_dense);
	return M_dense;
}

DenseVector<Real> 
RomSolverUtilities::reduceVector(DenseVector<Real> & V)
{
	DenseVector<Real> Vr;
	Q_dense.vector_mult_transpose(Vr, V);
	return Vr;
}

// ************************* Solution Update Functions ************************* //
DenseVector<Real> 
RomSolverUtilities::reconstructVector(DenseVector<Real> & Vr)
{
	DenseVector<Real> V;
	Q_dense.vector_mult(V, Vr);
	return V;
}

DenseVector<Real>
RomSolverUtilities::updateSolutionOnBoundary(DenseVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map)
{
	DenseVector<Real> U(S);
	for(std::map<int, std::pair<std::string, Real>>::iterator it = BC_map.begin(); it != BC_map.end(); ++it)
		U(it->first) = BC_map[it->first].second;
	return U;
}

DenseVector<Real>
RomSolverUtilities::updateSolutionOnBoundary(NumericVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map)
{
	for(std::map<int, std::pair<std::string, Real>>::iterator it = BC_map.begin(); it != BC_map.end(); ++it)
		S.set(it->first, BC_map[it->first].second);
	DenseVector<Real> U(S.size());
	for (int i=0; i<S.size(); i++)
		U(i) = S.el(i);
	return U;
}

void
RomSolverUtilities::updateSolutionOnBoundary_N(NumericVector<Real> & S, std::map<int, std::pair<std::string, Real>> BC_map)
{
	for(std::map<int, std::pair<std::string, Real>>::iterator it = BC_map.begin(); it != BC_map.end(); ++it)
		S.set(it->first, BC_map[it->first].second);
}

// ************************* Scaling Functions ************************* //

DenseVector<Real>
RomSolverUtilities::scaleVector(DenseVector<Real> V, DenseVector<Real> add_vec, DenseVector<Real> mult_vec, bool inverse, bool do_add, bool do_mult)
{
	DenseVector<Real> U(V.size());
	for (int i=0; i<V.size(); i++)
	{
		if (do_add and do_mult)
		{
			if (inverse)
			{
				if (mult_vec.el(i) != 0)
					U(i) = (V.el(i) - add_vec.el(i)) / mult_vec.el(i);
				else
					U(i) = V.el(i) - add_vec.el(i);
			}
			else
				U(i) = (V.el(i)*mult_vec.el(i)) + add_vec.el(i);
		}
		else if (do_add and !do_mult)
		{
			if (inverse)
			{
				U(i) = V.el(i) - add_vec.el(i);
			}
			else
				U(i) = V.el(i) + add_vec.el(i);
		}
		else if (!do_add and do_mult)
		{
			if (inverse)
			{
				if (mult_vec.el(i) != 0)
					U(i) = V.el(i) / mult_vec.el(i);
				else
					U(i) = 0.0;
			}
			else
				U(i) = V.el(i)*mult_vec.el(i);
		}
		else
			U(i) = V.el(i);
	}
	return U;
}


