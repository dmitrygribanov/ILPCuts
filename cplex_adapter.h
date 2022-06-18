#pragma once

#include "basic_utils.h"
#include <ilcplex/ilocplex.h>

namespace ILPCuts
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
		using Strategies = typename generator_base::Strategies;

		bool lazy_constraints_on = false;

		std::unique_ptr<generator_base> generator_ptr;

	private:
		std::vector<is_vec> cuts_Eigen;

		static constexpr auto to_IloInt = Utility::conv<int_type, IloInt>;
		static constexpr auto from_IloNum = Utility::conv<IloNum, int_type>;
		// Utility::ToSingle<int_type, IloInt> to_IloInt;

		std::tuple<is_mat, is_vec, rs_vec>
		transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars);

	public:
		void add_cuts(IloCplex cplex, IloNumVarArray vars);
		void generate_cuts(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars);
		void transform_cuts_to_ilo_constraints(IloCplex cplex, IloNumVarArray vars, IloRangeArray ilo_constr);
		void add_cuts_as_constraints(IloCplex cplex, IloModel model, IloRangeArray constrs, IloNumVarArray vars);
	};

	template <typename cutsGenT>
	std::tuple<typename CPlexAdapter<cutsGenT>::is_mat, typename CPlexAdapter<cutsGenT>::is_vec, typename CPlexAdapter<cutsGenT>::rs_vec>
	CPlexAdapter<cutsGenT>::transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
	{
		size_t n = vars.getSize();

		is_mat A_basis(n, n);
		rs_vec optLP(n);
		is_vec c(n);

		std::map<IloInt, size_t> var_index_map;
		for (size_t i = 0; i < n; ++i)
		{
			var_index_map[vars[i].getId()] = i;
		}

		{
			IloExpr::LinearIterator iter = obj.getLinearIterator();
			if (obj.getSense() == IloObjective::Sense::Maximize)
			{
				while (iter.ok())
				{
					c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum(iter.getCoef());
					++iter;
				}
			}
			else
			{
				// cout << "Min " << endl;
				while (iter.ok())
				{
					// cout << iter.getCoef() << endl;
					c.insertBack(var_index_map[iter.getVar().getId()]) = from_IloNum(-iter.getCoef());
					++iter;
				}
			}
		}

		{
			IloCplex::BasisStatusArray vars_basis(cplex.getEnv());
			IloCplex::BasisStatusArray constrs_basis(cplex.getEnv());
			cplex.getBasisStatuses(vars_basis, vars, constrs_basis, constrs);
			size_t basis_i = 0;

			using Trilpl = Triplet<int_type>;
			std::vector<Trilpl> triplets;

			for (size_t i = 0; i < constrs.getSize(); ++i)
			{
				IloRange constr = constrs[i];
				// cout << "LB, UB" << constr.getLB() << ", " << constr.getUB() << endl;
				/*if (constr.getLB() == constr.getUB()) {
					cout << "LB = UB" << endl;
					cin.ignore();
				}*/

				if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper)
				{
					IloExpr::LinearIterator iter = constr.getLinearIterator();
					while (iter.ok())
					{
						// cout << "AtUpper Row" << endl;
						triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], from_IloNum(iter.getCoef())));
						// A_basis(basis_i, var_index_map[iter.getVar().getId()]) = iter.getCoef();
						++iter;
					}
					++basis_i;
				}
				else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower)
				{
					IloExpr::LinearIterator iter = constr.getLinearIterator();
					while (iter.ok())
					{
						// cout << "AtLower Row" << endl;
						triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], from_IloNum(-iter.getCoef())));
						// A_basis(basis_i, var_index_map[iter.getVar().getId()]) = -iter.getCoef();
						++iter;
					}
					++basis_i;
				}
			}

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

		A_basis = generator_ptr->adjust_normal_cone(std::move(A_basis.transpose()), c);

		return std::make_tuple(A_basis, c, optLP);
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::add_cuts(IloCplex cplex, IloNumVarArray vars)
	{
		IloEnv env = cplex.getEnv();

		IloRangeArray ilo_cuts(env);

		transform_cuts_to_ilo_constraints(cplex, vars, ilo_cuts);

		cplex.addUserCuts(ilo_cuts);
		if (lazy_constraints_on)
			cplex.addLazyConstraints(ilo_cuts);

		ilo_cuts.end();
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::generate_cuts(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars)
	{
		is_mat A;
		rs_vec optLP;
		is_vec c;

		tie(A, c, optLP) = transform_to_corner_polyhedron(cplex, model, obj, constrs, vars);

		generator_ptr->direction_strategies.push_back(std::make_unique<typename Strategies::GreadyDirection>(c));

		cuts_Eigen = std::move(generator_ptr->generate_cuts(A, optLP));

		generator_ptr->direction_strategies.pop_back();
	}

	template <typename cutsGenT>
	void CPlexAdapter<cutsGenT>::transform_cuts_to_ilo_constraints(IloCplex cplex, IloNumVarArray vars, IloRangeArray ilo_constr)
	{
		IloEnv env = cplex.getEnv();

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
	void CPlexAdapter<cutsGenT>::add_cuts_as_constraints(IloCplex cplex, IloModel model, IloRangeArray constrs, IloNumVarArray vars)
	{
		IloEnv env = cplex.getEnv();

		IloRangeArray ilo_constraints(env);

		transform_cuts_to_ilo_constraints(cplex, vars, ilo_constraints);

		constrs.add(ilo_constraints);
		model.add(ilo_constraints);

		ilo_constraints.end();
	}
}
