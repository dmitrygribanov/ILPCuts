#pragma once

#include "base_module.h"
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplexx.h>
#include <map>

namespace ILPCuts
{
	string create_err_string(string prefix, CPXCENVptr env, int status, string suff = "")
	{
		stringstream res_msg;
		res_msg << prefix;

		char buff_C[CPXMESSAGEBUFSIZE];
		const char *msg_C = CPXXgeterrorstring(env, status, buff_C);

		if (msg_C)
			res_msg << msg_C;
		else
			res_msg << "CPlex Error " << status << ": unknown error code.\n";

		res_msg << suff;

		return res_msg.str();
	}

	map<string, IloCplex::IntParam> init_cplex_cut_type_names_map()
	{
		map<string, IloCplex::IntParam> res;

		res["clique"] = IloCplex::Param::MIP::Cuts::Cliques;
		res["cover"] = IloCplex::Param::MIP::Cuts::Covers;
		res["gomory"] = IloCplex::Param::MIP::Cuts::Gomory;
		res["disj"] = IloCplex::Param::MIP::Cuts::Disjunctive;
		res["flowc"] = IloCplex::Param::MIP::Cuts::FlowCovers;
		res["flowp"] = IloCplex::Param::MIP::Cuts::PathCut;
		res["gub"] = IloCplex::Param::MIP::Cuts::GUBCovers;
		res["liftp"] = IloCplex::Param::MIP::Cuts::LiftProj;
		res["mcf"] = IloCplex::Param::MIP::Cuts::MCFCut;
		res["mir"] = IloCplex::Param::MIP::Cuts::MIRCut;
		res["zeroh"] = IloCplex::Param::MIP::Cuts::ZeroHalfCut;
		res["impliedl"] = IloCplex::Param::MIP::Cuts::LocalImplied;
		res["impliedg"] = IloCplex::Param::MIP::Cuts::Implied;

		return res;
	}
	class FractionalElementException : public Exception
	{
	public:
		FractionalElementException()
			: Exception("ILPCuts lib error: fractional element in A-matrix has been found") {}

		FractionalElementException(std::string msg)
			: Exception(msg) {}
	};
	template <typename cutsGenT>
	class CPlexAdapter
	{
	public:
		using generator_base = cutsGenT;

		using int_type = typename generator_base::int_type;
		using rational_type = typename generator_base::rational_type;
		using is_mat = typename generator_base::is_mat;
		using rs_mat = typename generator_base::rs_mat;
		using is_vec = typename generator_base::is_vec;
		using rs_vec = typename generator_base::rs_vec;
		using Strategies = typename generator_base::Strategies;

		enum class IsLazy
		{
			on,
			off,
			hybrid
		};

		IsLazy is_lazy = IsLazy::off;

		std::unique_ptr<generator_base> generator_ptr;

	private:
		std::vector<is_vec> cuts_Eigen;

		static constexpr auto to_IloInt = Utility::conv<int_type, IloInt>;
		static constexpr auto from_IloNum_to_int = Utility::conv<IloNum, int_type>;
		static constexpr auto from_IloNum_to_rat = Utility::conv<IloNum, rational_type>;
		// static constexpr auto from_CPXINT_to_int = Utility::conv<CPXINT, rational_type>;

		std::tuple<is_mat, rs_vec, rs_vec>
		transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars);

		std::tuple<is_mat, rs_vec, rs_vec>
		transform_to_corner_polyhedron(CPXCENVptr env, CPXCLPptr lp);

		using Trilpl = Triplet<int_type>;
		void _fill_triplets_IloRow(vector<Trilpl> *triplets, IloExpr::LinearIterator *ilo_iter, int_type sign, size_t basis_i, const map<IloInt, size_t> &var_index_map, const IloRange &constr)
		{
			while (ilo_iter->ok())
			{
				rational_type coeff = from_IloNum_to_rat(ilo_iter->getCoef());

				int_type coeff_int = 0;

				try
				{
					coeff_int = generator_ptr->_rat_to_int(coeff);
				}
				catch (ToleranceErrorException)
				{
					stringstream msg;

					msg << "ILPCuts lib error: fracttional elemen in A has located;\n"
						<< "ConstrID = " << constr.getId() << "\n"
						<< "VarID = " << ilo_iter->getVar().getId() << "\n"
						<< "value = " << coeff;

					throw FractionalElementException(msg.str());
				}

				// rational_type rounded_coeff = round(coeff);
				// if (abs(coeff - rounded_coeff) >= generator_ptr->err_tolerance)
				// 	throw FractionalElementException();

				triplets->push_back(Trilpl(basis_i, var_index_map.at(ilo_iter->getVar().getId()), sign * coeff_int));

				++(*ilo_iter);
			}
			++basis_i;
		}

	public:
		logger &get_logger() { return generator_ptr->log; }
		logger::stream &out()
		{
			get_logger().out() << "ILPCuts-CPlex-Adapter: ";
			return get_logger().out();
		}
		void add_cuts(IloCplex cplex, IloNumVarArray vars);
		void add_cuts(CPXCENVptr env, void *cbdata, int wherefrom);
		void generate_cuts(is_mat A, rs_vec c, rs_vec optLP);
		void generate_cuts(CPXCENVptr env, CPXCLPptr lp)
		{
			is_mat A;
			rs_vec c, optLP;

			out() << "transforming data from ILO-C-Format to ILPCuts format ..." << endl;

			Timer timer;
			timer.start();

			tie(A, c, optLP) = transform_to_corner_polyhedron(env, lp);

			out() << "transformation time: " << timer.finish_str();

			generate_cuts(move(A), move(c), move(optLP));
		}
		void generate_cuts(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
		{
			is_mat A;
			rs_vec c, optLP;

			out() << "transforming data from ILO-C++-Format to ILPCuts format ..." << endl;

			Timer timer;
			timer.start();

			tie(A, c, optLP) = transform_to_corner_polyhedron(cplex, model, obj, constrs, vars);

			out() << "transformation time: " << timer.finish_str();

			generate_cuts(move(A), move(c), move(optLP));
		}
		void transform_cuts_to_ilo_constraints(IloEnv env, IloNumVarArray vars, IloRangeArray ilo_constr);

		void add_cuts_as_constraints(IloEnv env, IloNumVarArray vars, IloModel model, IloRangeArray constrs);
		void add_cuts_as_constraints(IloEnv env, IloNumVarArray vars, IloModel model);

		void clear_cuts()
		{
			cuts_Eigen.clear();
		}

		std::vector<is_vec> get_cuts() const { return cuts_Eigen; }
		void export_cuts(ofstream &out_file)
		{
			out_file << Utility::export_cuts(cuts_Eigen);
		}

		bool set_lazy(string str_par)
		{
			if (str_par == "on")
				is_lazy = IsLazy::on;
			else if (str_par == "off")
				is_lazy = IsLazy::off;
			else if (str_par == "hybrid")
				is_lazy == IsLazy::hybrid;
			else
				return false;

			return true;
		}
	};

	template <typename cutsGenT>
	std::tuple<typename CPlexAdapter<cutsGenT>::is_mat, typename CPlexAdapter<cutsGenT>::rs_vec, typename CPlexAdapter<cutsGenT>::rs_vec>
	CPlexAdapter<cutsGenT>::transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
	{
		size_t n = vars.getSize();

		is_mat A_basis(n, n);
		rs_vec optLP(n);
		rs_vec c(n);

		std::map<IloInt, size_t> var_index_map;
		for (size_t i = 0; i < n; ++i)
		{
			var_index_map[vars[i].getId()] = i;
		}

		out() << "reading objective..." << endl;

		{
			IloExpr::LinearIterator iter = obj.getLinearIterator();
			if (obj.getSense() == IloObjective::Sense::Maximize)
			{
				while (iter.ok())
				{
					c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum_to_rat(iter.getCoef());
					++iter;
				}
			}
			else
			{
				// cout << "Min " << endl;
				while (iter.ok())
				{
					// cout << iter.getCoef() << endl;
					c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum_to_rat(-iter.getCoef());
					++iter;
				}
			}
		}

		{
			out() << "reading active constraints..." << endl;

			IloCplex::BasisStatusArray vars_basis(cplex.getEnv());
			IloCplex::BasisStatusArray constrs_basis(cplex.getEnv());
			cplex.getBasisStatuses(vars_basis, vars, constrs_basis, constrs);
			size_t basis_i = 0;

			std::vector<Trilpl> triplets;

			for (size_t i = 0; i < constrs.getSize(); ++i)
			{
				IloRange constr = constrs[i];

				IloExpr::LinearIterator iter = constr.getLinearIterator();

				if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper)
				{
					_fill_triplets_IloRow(&triplets, &iter, 1, basis_i, var_index_map, constr);
					++basis_i;
				}
				else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower)
				{
					_fill_triplets_IloRow(&triplets, &iter, -1, basis_i, var_index_map, constr);
					++basis_i;
				}
			}

			out() << "reading active variables..." << endl;

			for (size_t i = 0; i < vars.getSize(); ++i)
			{
				IloNumVar var = vars[i];
				// cout << "LB, UB" << var.getLB() << ", " << var.getUB() << endl;
				if (vars_basis[i] == IloCplex::BasisStatus::AtUpper)
				{
					// cout << "AtUpper Var" << endl;
					triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], 1));
					// A_basis(basis_i, var_index_map[var.getId()]) = 1;
					++basis_i;
				}
				else if (vars_basis[i] == IloCplex::BasisStatus::AtLower)
				{
					// cout << "AtLower Var" << endl;
					triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], -1));
					// A_basis(basis_i, var_index_map[var.getId()]) = -1;
					++basis_i;
				}
			}

			// cin.ignore();

			A_basis.setFromTriplets(triplets.begin(), triplets.end());

			vars_basis.end();
			constrs_basis.end();
		}

		out() << "reading optimal LP solution (vertex)..." << endl;

		{
			IloNumArray cplex_optLP(cplex.getEnv());
			cplex.getValues(vars, cplex_optLP);

			for (size_t i = 0; i < n; ++i)
			{
				if (cplex_optLP[i] != 0)
					optLP.insertBack(i) = cplex_optLP[i];
			}

			cplex_optLP.end();
		}

		out() << "normal cone adjusting..." << endl;

		A_basis = generator_ptr->adjust_normal_cone(std::move(A_basis.transpose()), c);

		return std::make_tuple(A_basis, c, optLP);
	}

	template <typename cutsGenT>
	std::tuple<typename CPlexAdapter<cutsGenT>::is_mat, typename CPlexAdapter<cutsGenT>::rs_vec, typename CPlexAdapter<cutsGenT>::rs_vec>
	CPlexAdapter<cutsGenT>::transform_to_corner_polyhedron(CPXCENVptr env, CPXCLPptr lp)
	{
		int status = 0;

		CPXDIM n = 0, m = 0;
		n = CPXXgetnumcols(env, lp);
		m = CPXXgetnumrows(env, lp);

		is_mat A_basis(n, n);
		rs_vec optLP(n);
		rs_vec c(n);

		{
			out() << "reading objective..." << endl;

			int sence = CPXXgetobjsen(env, lp);
			if (sence != CPX_MIN && sence != CPX_MAX)
			{
				throw Exception("ILPCuts lib error: the LP problem's optimizaton sence var is = 0, so it is not defined...\n");
			}
			// CPX_MAX == -1 and CPX_MIN == 1 : we need to change into opposite values
			sence *= -1;

			vector<double> C_obj(n, 0.0);
			status = CPXXgetobj(env, lp, C_obj.data(), 0, n - 1);
			if (status)
			{
				throw Exception(create_err_string("ILPCuts lib error: failed to read objective using CPXXgetobj...\n",
												  env, status));
			}

			for (size_t i = 0; i < C_obj.size(); ++i)
			{
				if (abs(C_obj[i]) > numeric_limits<double>::epsilon())
				{
					c.insertBack(i) = from_IloNum_to_rat(sence * C_obj[i]);
				}
			}
		}

		{
			out() << "reading basis information..." << endl;

			vector<int> cols_status(n, 0);
			vector<int> rows_status(m, 0);
			std::vector<Trilpl> triplets;
			size_t basis_i = 0;

			status = CPXXgetbase(env, lp, cols_status.data(), rows_status.data());
			if (status)
			{
				throw Exception(create_err_string("ILPCuts lib error: failed to read base using CPXXgetbase...\n",
												  env, status));
			}

			out() << "reading active constraints..." << endl;
			{
				CPXNNZ nonzeros_C = 0, surplus_C = 0;
				status = CPXXgetrows(env, lp, &nonzeros_C, nullptr, nullptr, nullptr, 0, &surplus_C, 0, m - 1);
				if (status != CPXERR_NEGATIVE_SURPLUS)
				{
					stringstream msg_suff;
					msg_suff << "surplus is " << surplus_C << "; nonzeros is " << nonzeros_C << endl;
					throw Exception(create_err_string("ILPCuts lib error: failed to read nonzeros number in rows...\n",
													  env, status, msg_suff.str()));
				}
				nonzeros_C = -surplus_C;

				vector<CPXNNZ> nonzeros_beg_C(m, 0);
				vector<CPXDIM> mat_idxcs_C(nonzeros_C, 0);
				vector<double> mat_vals_C(nonzeros_C, 0.0);
				status = CPXXgetrows(env, lp, &nonzeros_C, nonzeros_beg_C.data(),
									 mat_idxcs_C.data(), mat_vals_C.data(), nonzeros_C,
									 &surplus_C, 0, m - 1);
				if (status && surplus_C != 0)
				{
					stringstream msg_suff;
					msg_suff << "surplus is " << surplus_C << "; nonzeros is " << nonzeros_C << endl;
					throw Exception(create_err_string("ILPCuts lib error: failed to read nonzeros from CPLex rows...\n",
													  env, status, msg_suff.str()));
				}

				// size_t bc = 0;
				// for (size_t base_row_idx = 0; base_row_idx < rows_status.size(); ++base_row_idx)
				// {

				// 	if (rows_status[base_row_idx] == CPX_AT_UPPER)
				// 		++bc;
				// 	else if (rows_status[base_row_idx] == CPX_AT_LOWER)
				// 		++bc;
				// }

				// out() << "BC = " << bc << endl;

				for (size_t base_row_idx = 0; base_row_idx < rows_status.size(); ++base_row_idx)
				{
					int mult = 0;
					if (rows_status[base_row_idx] == CPX_AT_UPPER)
						mult = 1;
					else if (rows_status[base_row_idx] == CPX_AT_LOWER)
						mult = -1;

					if (mult != 0)
					{
						CPXNNZ row_st_C = 0, row_end_C = 0;
						row_st_C = nonzeros_beg_C[base_row_idx];
						if (base_row_idx == rows_status.size() - 1)
							row_end_C = nonzeros_C - 1;
						else
							row_end_C = nonzeros_beg_C[base_row_idx + 1];

						for (CPXNNZ nonzeros_idx = row_st_C; nonzeros_idx < row_end_C; ++nonzeros_idx)
						{
							CPXDIM col_idx = mat_idxcs_C[nonzeros_idx];
							double val = mult * mat_vals_C[nonzeros_idx];
							int_type val_int = 0;

							// out() << "(" << basis_i << ", " << col_idx << ", " << val << ")" << endl;

							try
							{
								val_int = generator_ptr->_rat_to_int(val);
							}
							catch (ToleranceErrorException)
							{
								stringstream msg;

								msg << "ILPCuts lib error: fracttional elemen in A has located;\n"
									<< "ConstrNum = " << base_row_idx << "\n"
									<< "VarNum = " << col_idx << "\n"
									<< "value = " << val;

								throw FractionalElementException(msg.str());
							}

							triplets.emplace_back(Trilpl(basis_i, col_idx, val_int));
						}

						++basis_i;
					}
				}
			}

			out() << "reading active variables..." << endl;

			for (size_t base_var_idx = 0; base_var_idx < cols_status.size(); ++base_var_idx)
			{
				int sign = 0;
				if (cols_status[base_var_idx] == CPX_AT_UPPER)
					sign = 1;
				else if (cols_status[base_var_idx] == CPX_AT_LOWER)
					sign = -1;

				if (sign != 0)
				{
					triplets.push_back(Trilpl(basis_i, base_var_idx, sign));
					++basis_i;
				}
			}

			A_basis.setFromTriplets(triplets.begin(), triplets.end());

			// cout << "basis_i = " << basis_i << endl;
			// cin.ignore();
		}

		{
			out() << "reading optimal LP solution (vertex)..." << endl;

			vector<double> x_C(n, 0.0);
			status = CPXXsolution(env, lp, nullptr, nullptr, x_C.data(), nullptr, nullptr, nullptr);
			if (status)
			{
				throw Exception(create_err_string("ILPCuts lib error: failed to read lp-solution...\n",
												  env, status));
			}
			for (size_t i = 0; i < n; ++i)
			{
				if (abs(x_C[i]) > numeric_limits<double>::epsilon())
				{
					optLP.insertBack(i) = from_IloNum_to_rat(x_C[i]);
				}
			}
		}

		out() << "normal cone adjusting..." << endl;

		// cout << A_basis << endl;
		// cout << A_basis.cols() << ",  " << A_basis.rows() << endl;

		A_basis = generator_ptr->adjust_normal_cone(std::move(A_basis.transpose()), c);

		return std::make_tuple(A_basis, c, optLP);
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::add_cuts(CPXCENVptr env, void *cbdata, int wherefrom)
	{
		int status = 0;

		for (size_t i = 0; i < cuts_Eigen.size(); ++i)
		{
			const auto &cut = cuts_Eigen[i];
			CPXNNZ nonzeros_C = 0;
			int_type rhs_C = 0;

			nonzeros_C = cut.nonZeros();
			rhs_C = cut.coeff(cut.size() - 1);
			if (rhs_C != 0)
				--nonzeros_C;

			vector<CPXDIM> nonzero_idxs_C;
			vector<double> nonzero_vals_C;
			for (typename SparseVector<int_type>::InnerIterator it(cut); it; ++it)
			{
				if (it.index() != cut.size() - 1)
				{
					IloInt int_val = to_IloInt(it.value());
					nonzero_idxs_C.push_back(it.index());
					nonzero_vals_C.push_back(int_val);
					// myData->log->out() << "VALLL = " << it.value() << endl;
				}
			}

			// for (size_t i = 0; i < nonzeros_C; i++)
			// {
			//     cout << nonzero_vals_C[i] << " ";
			// }
			// cout << endl;

			status = CPXXcutcallbackaddlocal(env, cbdata, wherefrom, nonzeros_C, rhs_C, 'L',
											 nonzero_idxs_C.data(), nonzero_vals_C.data());
			// FOR GLOBAL CUTS: could be CPX_USECUT_FORCE, CPX_USECUT_PURGE or CPX_USECUT_FILTER
			if (status)
			{
				stringstream msg_suff;
				msg_suff << "the cut sequential number is " << i << endl;
				throw(Exception(create_err_string("fail to add an user cut...\n", env, status,
												  msg_suff.str())));
			}
		}
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::add_cuts(IloCplex cplex, IloNumVarArray vars)
	{
		IloEnv env = cplex.getEnv();

		IloRangeArray ilo_cuts(env);

		transform_cuts_to_ilo_constraints(env, vars, ilo_cuts);

		if (is_lazy == IsLazy::off)
		{
			cplex.addUserCuts(ilo_cuts);
		}
		else if (is_lazy == IsLazy::on)
		{
			cplex.addLazyConstraints(ilo_cuts);
		}
		else
		{
			cplex.addUserCuts(ilo_cuts);
			cplex.addLazyConstraints(ilo_cuts);
		}

		ilo_cuts.end();
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::generate_cuts(is_mat A, rs_vec c, rs_vec optLP)
	{
		generator_ptr->direction_strategies.push_back(std::make_unique<typename Strategies::GreadyDirection>(c));

		Timer timer;
		timer.start();
		out() << "cutting plane generation is started..." << endl;

		cuts_Eigen = generator_ptr->generate_cuts(A, optLP);

		out() << "cutting plane generation is finished,\n"
			  << "total time: " << timer.finish_str() << "\n"
			  << "number of added cuts: " << cuts_Eigen.size() << "\n"
			  << endl;

		generator_ptr->direction_strategies.pop_back();
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::transform_cuts_to_ilo_constraints(IloEnv env, IloNumVarArray vars, IloRangeArray ilo_constr)
	{
		size_t n = vars.getSize();

		for (size_t i = 0; i < cuts_Eigen.size(); ++i)
		{
			IloExpr lhs = IloExpr(env);
			IloNum rhs = 0;

			for (typename is_vec::InnerIterator it(cuts_Eigen[i]); it; ++it)
			{
				if (it.index() < n)
					lhs += vars[it.index()] * to_IloInt(it.value());
				else
					rhs = to_IloInt(it.value());
			}

			ilo_constr.add(IloRange(env, lhs, rhs));
		}
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::add_cuts_as_constraints(IloEnv env, IloNumVarArray vars, IloModel model, IloRangeArray constrs)
	{
		IloRangeArray ilo_constraints(env);

		transform_cuts_to_ilo_constraints(env, vars, ilo_constraints);

		constrs.add(ilo_constraints);
		model.add(ilo_constraints);

		ilo_constraints.end();
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::add_cuts_as_constraints(IloEnv env, IloNumVarArray vars, IloModel model)
	{
		IloRangeArray ilo_constraints(env);

		transform_cuts_to_ilo_constraints(env, vars, ilo_constraints);

		model.add(ilo_constraints);

		ilo_constraints.end();
	}
}
