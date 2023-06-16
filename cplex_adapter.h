#pragma once

#include "cplex_utils.h"
#include "base_algo.h"
namespace ILPCuts
{
	namespace ModuleCPlex
	{
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
			using Strategies = typename generator_base::ConeCoverer::Strategies;
			using cut_info = typename generator_base::cut_info;

			enum class IsLazy
			{
				on,
				off,
				hybrid
			};

			IsLazy is_lazy = IsLazy::off;

			std::unique_ptr<generator_base> generator_ptr;

		private:
			static constexpr auto to_IloInt = conv<int_type, IloInt>::make;
			static constexpr auto from_IloNum_to_int = conv<IloNum, int_type>::make;
			static constexpr auto from_IloNum_to_rat = conv<IloNum, rational_type>::make;
			// static constexpr auto from_CPXINT_to_int = conv<CPXINT, rational_type>;

			map<IloInt, size_t> init_var_index_map(IloNumVarArray vars) const
			{
				size_t n = vars.getSize();

				map<IloInt, size_t> res;
				for (size_t i = 0; i < n; ++i)
				{
					res[vars[i].getId()] = i;
				}
				return res;
			}

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

			std::tuple<is_mat, rs_vec, rs_vec>
			extract_problem_info(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars);

			std::tuple<is_mat, rs_vec, rs_vec>
			extract_problem_info(CPXCENVptr env, CPXCLPptr lp);

			void generate_cuts(CPXCENVptr env, CPXCLPptr lp)
			{
				is_mat A;
				rs_vec c, optLP;

				out() << "transforming data from ILO-C-Format to ILPCuts format ..." << endl;

				Timer timer;
				timer.start();

				tie(A, c, optLP) = extract_problem_info(env, lp);

				out() << "transformation time: " << timer.finish_str();

				generate_cuts(move(A), move(c), move(optLP));
			}
			void generate_cuts(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
			{
				is_mat A;
				rs_vec c, optLP;

				out() << "transforming data from ILO-C++-Format to ILPCuts format ..." << endl;

				Timer timer;
				timer.start();

				tie(A, c, optLP) = extract_problem_info(cplex, obj, constrs, vars);

				out() << "transformation time: " << timer.finish_str();

				generate_cuts(move(A), move(c), move(optLP));
			}
			void transform_cuts_to_ilo_constraints(IloEnv env, IloNumVarArray vars, IloRangeArray ilo_constr);

			IloRangeArray create_IloConstraints_from_cuts(IloEnv env, IloNumVarArray vars)
			{
				IloRangeArray ilo_constraints(env);

				transform_cuts_to_ilo_constraints(env, vars, ilo_constraints);

				return ilo_constraints;
			}

			void clear_cuts()
			{
				generator_ptr->clear_cuts();
			}

			const std::vector<cut_info> &get_cuts() const { return generator_ptr->get_cuts(); }

			void export_cuts(ofstream &out_file)
			{
				generator_ptr->export_cuts(out_file);
			}

			void import_cuts(ifstream &in_file)
			{
				generator_ptr->import_cuts(in_file);
			}

			void check_cuts_validity()
			{
				generator_ptr->check_cuts_validity();
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

			IloObjective apply_objective(const is_vec &new_obj, IloModel model, IloNumVarArray vars)
			{
				IloEnv env = model.getEnv();

				IloExpr expr(env);
				for (typename is_vec::InnerIterator it(new_obj); it; ++it)
				{
					expr += vars[it.index()] * to_IloInt(it.value());
				}

				IloObjective res_obj(env, expr, IloObjective::Sense::Maximize);

				model.add(res_obj);
				return res_obj;
			}
		};

		template <typename cutsGenT>
		std::tuple<typename CPlexAdapter<cutsGenT>::is_mat, typename CPlexAdapter<cutsGenT>::rs_vec, typename CPlexAdapter<cutsGenT>::rs_vec>
		CPlexAdapter<cutsGenT>::extract_problem_info(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
		{
			size_t n = vars.getSize();

			is_mat A_basis(n, n);
			rs_vec optLP(n);
			rs_vec c(n);

			auto var_index_map = init_var_index_map(vars);

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

					if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper ||
						constrs_basis[i] == IloCplex::BasisStatus::AtLower && constr.getLB() == -numeric_limits<IloNum>::infinity())
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
				if (basis_i != n)
				{
					stringstream msg;
					msg << "ILPCuts error: failed to read IloCplex problem\n"
						<< "size of the basis (" << basis_i << "is not equal to n (" << n << ")";
					throw(Exception(msg.str()));
				}

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

			A_basis = generator_ptr->adjust_normal_cone(std::move(A_basis.transpose()), c).transpose().eval();

			return std::make_tuple(A_basis, c, optLP);
		}

		template <typename cutsGenT>
		std::tuple<typename CPlexAdapter<cutsGenT>::is_mat, typename CPlexAdapter<cutsGenT>::rs_vec, typename CPlexAdapter<cutsGenT>::rs_vec>
		CPlexAdapter<cutsGenT>::extract_problem_info(CPXCENVptr env, CPXCLPptr lp)
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

			A_basis = generator_ptr->adjust_normal_cone(std::move(A_basis.transpose()), c).transpose().eval();

			return std::make_tuple(A_basis, c, optLP);
		}

		template <typename cutsGenT>
		void CPlexAdapter<cutsGenT>::add_cuts(CPXCENVptr env, void *cbdata, int wherefrom)
		{
			int status = 0;

			const auto &cuts_Eigen = get_cuts();

			for (size_t i = 0; i < cuts_Eigen.size(); ++i)
			{
				const auto &cut = cuts_Eigen[i];
				CPXNNZ nonzeros_C = 0;

				nonzeros_C = cut.ray.nonZeros();

				vector<CPXDIM> nonzero_idxs_C;
				vector<double> nonzero_vals_C;
				for (typename SparseVector<int_type>::InnerIterator it(cut.ray); it; ++it)
				{
					IloInt int_val = to_IloInt(it.value());
					nonzero_idxs_C.push_back(it.index());
					nonzero_vals_C.push_back(int_val);
				}

				// for (size_t i = 0; i < nonzeros_C; i++)
				// {
				//     cout << nonzero_vals_C[i] << " ";
				// }
				// cout << endl;

				status = CPXXcutcallbackaddlocal(env, cbdata, wherefrom, nonzeros_C, cut.rhs, 'L',
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
			generator_ptr->cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename Strategies::GreadyDirection>(c));

			Timer timer;
			timer.start();
			out() << "cutting plane generation is started..." << endl;

			generator_ptr->LocalInstance.A = std::move(A);
			generator_ptr->LocalInstance.c = std::move(c);
			generator_ptr->LocalInstance.v = std::move(optLP);
			generator_ptr->generate_cuts();

			out() << "cutting plane generation is finished,\n"
				  << "total time: " << timer.finish_str() << "\n"
				  << "number of added cuts: " << get_cuts().size() << "\n"
				  << endl;

			generator_ptr->cone_coverer_ptr->direction_strategies.pop_back();
		}

		template <typename cutsGenT>
		void CPlexAdapter<cutsGenT>::transform_cuts_to_ilo_constraints(IloEnv env, IloNumVarArray vars, IloRangeArray ilo_constr)
		{
			size_t n = vars.getSize();

			const auto &cuts_Eigen = get_cuts();

			for (size_t i = 0; i < cuts_Eigen.size(); ++i)
			{
				IloExpr lhs = IloExpr(env);
				IloNum rhs = 0;

				for (typename is_vec::InnerIterator it(cuts_Eigen[i].ray); it; ++it)
				{
					lhs += vars[it.index()] * to_IloInt(it.value());
				}

				ilo_constr.add(IloRange(env, lhs, cuts_Eigen[i].rhs));
			}
		}

		template <typename cutsGenT>
		bool generate_cuts_with_time_limit(CPlexAdapter<cutsGenT> &cuts_gen_adapter, double time_limit, IloCplex solver, IloObjective obj, IloRangeArray rngs, IloNumVarArray vars)
		{
			using CPlexAdapter = CPlexAdapter<cutsGenT>;

			if (time_limit > 0.0)
			{
				// chose correct overloaded member "generate_cuts" from "CutsGenAdapter"
				auto ptr_generate_cuts = static_cast<void (CPlexAdapter::*)(IloCplex, IloObjective, IloRangeArray, IloNumVarArray)>(&CPlexAdapter::generate_cuts);

				// bind it to a compact function of the type "void()"
				auto bind_generate_cuts = bind(ptr_generate_cuts, &cuts_gen_adapter, solver, obj, rngs, vars);

				return run_async_function_with_time_limit(time_limit, bind_generate_cuts);
			}
			else
			{
				cuts_gen_adapter.generate_cuts(solver, obj, rngs, vars);

				return true;
			}
		}
	} // namespace ModuleCPlex
} // namespace ILPCuts
