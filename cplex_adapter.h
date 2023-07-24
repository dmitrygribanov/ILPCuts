#pragma once

#include "cuts_gen_base.h"
#include "cplex_utils.h"
namespace ILPCuts
{
	namespace ModuleCPlex
	{
		class CPlexAdapter
		{
		public:
			using generator_base = cuts_gen_base;

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

			using Trilpl = Triplet<rational_type>;
			void _fill_triplets_IloRow(vector<Trilpl> *triplets, IloExpr::LinearIterator *ilo_iter, rational_type sign, size_t basis_i, const map<IloInt, size_t> &var_index_map, const IloRange &constr)
			{
				while (ilo_iter->ok())
				{
					rational_type coeff = from_IloNum_to_rat(ilo_iter->getCoef());

					if (ilo_iter->getVar().getType() != IloNumVar::Float)
					{
						if (is_fract(coeff, generator_ptr->err_tolerance))
						{
							throw FractionalElementException(constr.getId(), ilo_iter->getVar().getId(), coeff);
						}
					}

					triplets->push_back(Trilpl(basis_i, var_index_map.at(ilo_iter->getVar().getId()), coeff * sign));

					// int_type coeff_int = 0;

					// try
					// {
					// 	coeff_int = generator_ptr->rat_to_int(coeff);
					// }
					// catch (ToleranceErrorException)
					// {
					// 	throw FractionalElementException(constr.getId(), ilo_iter->getVar().getId(), coeff);
					// }

					// triplets->push_back(Trilpl(basis_i, var_index_map.at(ilo_iter->getVar().getId()), sign * coeff_int));

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
			void generate_cuts(local_instance inst);

			local_instance extract_problem_info(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars);

			local_instance extract_problem_info(CPXCENVptr env, CPXCLPptr lp);

			void generate_cuts(CPXCENVptr env, CPXCLPptr lp)
			{
				is_mat A;
				rs_vec c, optLP;

				out() << "transforming data from ILO-C-Format to ILPCuts format ..." << endl;

				Timer timer;
				timer.start();

				local_instance inst = extract_problem_info(env, lp);

				out() << "transformation time: " << timer.finish_str();

				generate_cuts(inst);
			}
			void generate_cuts(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
			{
				is_mat A;
				rs_vec c, optLP;

				out() << "transforming data from ILO-C++-Format to ILPCuts format ..." << endl;

				Timer timer;
				timer.start();

				local_instance inst = extract_problem_info(cplex, obj, constrs, vars);

				out() << "transformation time: " << timer.finish_str();

				generate_cuts(inst);
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

			void check_pure_int_cuts_validity()
			{
				generator_ptr->check_pure_int_cuts_validity();
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

		local_instance CPlexAdapter::extract_problem_info(IloCplex cplex, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
		{
			size_t n = vars.getSize();
			local_instance resInst(n);

			auto var_index_map = init_var_index_map(vars);

			out() << "reading var.-types..." << endl;
			for (size_t i = 0; i < vars.getSize(); ++i)
			{
				IloNumVar var = vars[i];

				if (var.getType() == IloNumVar::Float)
				{
					resInst.ILP_var_types[var_index_map[var.getId()]] = local_instance::var_ILP_type::Rat;
				}
				else
				{
					resInst.ILP_var_types[var_index_map[var.getId()]] = local_instance::var_ILP_type::Int;
				}
				// cout << "LB, UB" << var.getLB() << ", " << var.getUB() << endl;
			}

			out() << "reading objective..." << endl;

			{
				IloExpr::LinearIterator iter = obj.getLinearIterator();
				if (obj.getSense() == IloObjective::Sense::Maximize)
				{
					while (iter.ok())
					{
						resInst.c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum_to_rat(iter.getCoef());
						++iter;
					}
				}
				else
				{
					while (iter.ok())
					{
						resInst.c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum_to_rat(-iter.getCoef());
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
						_fill_triplets_IloRow(&triplets, &iter, 1.0, basis_i, var_index_map, constr);
						++basis_i;
					}
					else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower)
					{
						_fill_triplets_IloRow(&triplets, &iter, -1.0, basis_i, var_index_map, constr);
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
						triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], 1.0));
						// A_basis(basis_i, var_index_map[var.getId()]) = 1;
						++basis_i;
					}
					else if (vars_basis[i] == IloCplex::BasisStatus::AtLower)
					{
						// cout << "AtLower Var" << endl;
						triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], -1.0));
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

				resInst.A.setFromTriplets(triplets.begin(), triplets.end());

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
						resInst.v.insertBack(i) = cplex_optLP[i];
				}

				cplex_optLP.end();
			}

			{
				out() << "reading LB and UB for variables..." << endl;

				out() << "reading var.-types..." << endl;
				for (size_t i = 0; i < vars.getSize(); ++i)
				{
					IloNumVar var = vars[i];

					size_t var_idxs = var_index_map[var.getId()];
					resInst.varsLB.insertBack(var_idxs) = from_IloNum_to_rat(var.getLb());
					resInst.varsUB.insertBack(var_idxs) = from_IloNum_to_rat(var.getUb());
				}
			}

			out() << "normal cone adjusting..." << endl;

			resInst.A = generator_ptr->adjust_normal_cone(std::move(resInst.A.transpose()), resInst.c).transpose().eval();

			return resInst;
		}

		local_instance CPlexAdapter::extract_problem_info(CPXCENVptr env, CPXCLPptr lp)
		{
			int status = 0;

			CPXDIM n = 0, m = 0;
			n = CPXXgetnumcols(env, lp);
			m = CPXXgetnumrows(env, lp);

			local_instance resInst(n);

			{
				out() << "reading var.-types..." << endl;
				std::vector<char> var_types_C(n, 0);
				int status = CPXXgetctype(env, lp, var_types_C.data(), 0, n);
				if (status)
				{
					throw Exception(create_err_string("ILPCuts lib error: failed to read var.-types using CPXXgetctype...\n",
													  env, status));
				}
				for (size_t i = 0; i < n; ++i)
				{
					if (var_types_C[i] == CPX_CONTINUOUS)
					{
						resInst.ILP_var_types[i] = local_instance::var_ILP_type::Rat;
					}
					else if (var_types_C[i] == CPX_BINARY || var_types_C[i] == CPX_INTEGER)
					{
						resInst.ILP_var_types[i] = local_instance::var_ILP_type::Int;
					}
					else
					{
						throw SemiVariableException(i, var_types_C[i]);
					}
				}
			}

			{
				out() << "reading LB and UB for variables..." << endl;

				std::vector<double> LBses_C(n), UBses_C(n);

				int status = CPXXgetlb(env, lp, LBses_C.data(), 0, n);
				if (status)
				{
					throw Exception(create_err_string("ILPCuts lib error: failed to read var's LB'es using CPXXgetlb...\n",
													  env, status));
				}

				status = CPXXgetub(env, lp, UBses_C.data(), 0, n);
				if (status)
				{
					throw Exception(create_err_string("ILPCuts lib error: failed to read var's UB'es using CPXXgetub...\n",
													  env, status));
				}

				for (size_t i = 0; i < n; ++i)
				{
					resInst.varsLB.insertBack(i) = from_IloNum_to_rat(LBses_C[i]);
					resInst.varsUB.insertBack(i) = from_IloNum_to_rat(UBses_C[i]);
				}
			}

			{
				out() << "reading objective..." << endl;

				int sence = CPXXgetobjsen(env, lp);
				if (sence != CPX_MIN && sence != CPX_MAX)
				{
					throw Exception("ILPCuts lib error: the LP problem's optimizaton sence var is = 0, so it is not defined...\n");
				}

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
						resInst.c.insertBack(i) = from_IloNum_to_rat(sence * C_obj[i]);
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
						double mult = 0;
						if (rows_status[base_row_idx] == CPX_AT_UPPER)
							mult = 1.0;
						else if (rows_status[base_row_idx] == CPX_AT_LOWER)
							mult = -1.0;

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

								if (resInst.ILP_var_types[col_idx] == local_instance::var_ILP_type::Int)
								{
									if (is_fract(val, generator_ptr->err_tolerance))
									{
										throw FractionalElementException(base_row_idx, col_idx, val);
									}
								}

								// int_type val_int = 0;

								// try
								// {
								// 	val_int = generator_ptr->rat_to_int(val);
								// }
								// catch (ToleranceErrorException)
								// {
								// 	throw FractionalElementException(base_row_idx, col_idx, val);
								// }

								triplets.emplace_back(Trilpl(basis_i, col_idx, val));
							}

							++basis_i;
						}
					}
				}

				out() << "reading active variables..." << endl;

				for (size_t base_var_idx = 0; base_var_idx < cols_status.size(); ++base_var_idx)
				{
					rational_type sign = 0.0;
					if (cols_status[base_var_idx] == CPX_AT_UPPER)
						sign = 1.0;
					else if (cols_status[base_var_idx] == CPX_AT_LOWER)
						sign = -1.0;

					if (sign != 0.0)
					{
						triplets.push_back(Trilpl(basis_i, base_var_idx, sign));
						++basis_i;
					}
				}

				resInst.A.setFromTriplets(triplets.begin(), triplets.end());
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
						resInst.v.insertBack(i) = from_IloNum_to_rat(x_C[i]);
					}
				}
			}

			out() << "normal cone adjusting..." << endl;

			resInst.A = generator_ptr->adjust_normal_cone(std::move(resInst.A.transpose()), resInst.c).transpose().eval();

			return resInst;
		}

		void CPlexAdapter::add_cuts(CPXCENVptr env, void *cbdata, int wherefrom)
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
				for (typename SparseVector<rational_type>::InnerIterator it(cut.ray); it; ++it)
				{
					IloInt int_val = to_IloInt(it.value());
					nonzero_idxs_C.push_back(it.index());
					nonzero_vals_C.push_back(int_val);
				}

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

		void CPlexAdapter::add_cuts(IloCplex cplex, IloNumVarArray vars)
		{
			IloEnv env = cplex.getEnv();

			if (get_cuts().size() != 0)
			{

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
		}

		void CPlexAdapter::generate_cuts(local_instance inst)
		{
			Timer timer;
			timer.start();
			out() << "cutting plane generation is started..." << endl;

			generator_ptr->LocalInstance = std::move(inst);
			generator_ptr->generate_cuts();

			out() << "cutting plane generation is finished,\n"
				  << "total time: " << timer.finish_str() << "\n"
				  << "number of added cuts: " << get_cuts().size() << "\n"
				  << endl;
		}

		void CPlexAdapter::transform_cuts_to_ilo_constraints(IloEnv env, IloNumVarArray vars, IloRangeArray ilo_constr)
		{
			size_t n = vars.getSize();

			const auto &cuts_Eigen = get_cuts();

			for (size_t i = 0; i < cuts_Eigen.size(); ++i)
			{
				IloExpr lhs = IloExpr(env);
				IloNum rhs = 0;

				for (typename rs_vec::InnerIterator it(cuts_Eigen[i].ray); it; ++it)
				{
					lhs += vars[it.index()] * it.value();
				}

				ilo_constr.add(IloRange(env, lhs, cuts_Eigen[i].rhs));
			}
		}

		bool generate_cuts_with_time_limit(CPlexAdapter &cuts_gen_adapter, double time_limit, IloCplex solver, IloObjective obj, IloRangeArray rngs, IloNumVarArray vars)
		{
			using CPlexAdapter = CPlexAdapter;

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
