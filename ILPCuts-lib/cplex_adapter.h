#pragma once

#include"base_module.h"

namespace ILPCuts {

    template<typename cutsGeneratorT>
    class CPlexAdapter :public cutsGeneratorT {
    public:

        using base = cutsGeneratorT;

        using typename base::int_type;
        using typename base::rational_type;
        using typename base::is_mat;
        using typename base::rs_mat;
        using typename base::is_vec;
        using typename base::rs_vec;

        using base::tolerance;
        using base::_rat_to_int;
        using base::adjust_normal_cone;

        std::vector<is_vec> cuts_Eigen;

        CPlexAdapter(std::ostream* outs, std::ostream* errs)
            :base(outs, errs) {}

        CPlexAdapter() :base() {}

        void add_cuts(IloCplex cplex, IloNumVarArray vars) {
            IloEnv env = cplex.getEnv();

            IloConstraintArray ilo_cuts(env);

            size_t n = vars.getSize();

            for (size_t i = 0; i < cuts_Eigen.size(); ++i) {
                IloExpr lhs = IloExpr(env);
                IloNum rhs = 0;

                for (typename is_vec::InnerIterator it(cuts_Eigen[i]); it; ++it) {
                    if (it.index() < n)
                        lhs += vars[it.index()] * it.value().get_si();
                    else
                        rhs = it.value().get_si();
                }

                ilo_cuts.add(IloRange(env, lhs, rhs));
            }

            cplex.addUserCuts(ilo_cuts);
            //cplex.addLazyConstraints(ilo_cuts);

            ilo_cuts.end();
        }

        void generate_cuts(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars) {
            is_mat A;
            rs_vec optLP;
            is_vec c;

            tie(A, c, optLP) = transform_to_corner_polyhedron(cplex, model, obj, constrs, vars);

            base::direction_strategies.push_back(std::make_unique<typename base::GreadyDirection>(c));

            cuts_Eigen = std::move(base::generate_cuts(A, optLP));
        }

        std::tuple<is_mat, is_vec, rs_vec>
            transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars) {
            size_t n = vars.getSize();

            typename base::is_mat A_basis(n, n);
            rs_vec optLP(n), c(n);

            std::map<IloInt, size_t> var_index_map;
            for (size_t i = 0; i < n; ++i) {
                var_index_map[vars[i].getId()] = i;
            }

            {
                IloExpr::LinearIterator iter = obj.getLinearIterator();
                if (obj.getSense() == IloObjective::Sense::Maximize)
                {
                    //cout << "Max " << endl;
                    while (iter.ok()) {
                        //cout << iter.getCoef() << endl;
                        c.insertBack(var_index_map[iter.getVar().getId()]) = iter.getCoef();
                        ++iter;
                    }
                }
                else {
                    //cout << "Min " << endl; 
                    while (iter.ok()) {
                        //cout << iter.getCoef() << endl;
                        c.insertBack(var_index_map[iter.getVar().getId()]) = -iter.getCoef();
                        ++iter;
                    }
                }

                /*cout << c << endl;
                cin.ignore();
                */
            }

            {
                IloCplex::BasisStatusArray vars_basis(cplex.getEnv());
                IloCplex::BasisStatusArray constrs_basis(cplex.getEnv());
                cplex.getBasisStatuses(vars_basis, vars, constrs_basis, constrs);
                size_t basis_i = 0;

                using Trilpl = Triplet<rational_type>;
                std::vector<Trilpl> triplets;

                for (size_t i = 0; i < constrs.getSize(); ++i) {
                    IloRange constr = constrs[i];
                    //cout << "LB, UB" << constr.getLB() << ", " << constr.getUB() << endl;
                    /*if (constr.getLB() == constr.getUB()) {
                        cout << "LB = UB" << endl;
                        cin.ignore();
                    }*/

                    if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            //cout << "AtUpper Row" << endl;
                            triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], iter.getCoef()));
                            //A_basis(basis_i, var_index_map[iter.getVar().getId()]) = iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                    else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            //cout << "AtLower Row" << endl;
                            triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], -iter.getCoef()));
                            //A_basis(basis_i, var_index_map[iter.getVar().getId()]) = -iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                }
                //cin.ignore();

                for (size_t i = 0; i < vars.getSize(); ++i) {
                    IloNumVar var = vars[i];
                    //cout << "LB, UB" << var.getLB() << ", " << var.getUB() << endl;
                    if (vars_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        //cout << "AtUpper Var" << endl;
                        triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], -1));
                        //A_basis(basis_i, var_index_map[var.getId()]) = 1;
                        ++basis_i;
                    }
                    else if (vars_basis[i] == IloCplex::BasisStatus::AtLower) {
                        //cout << "AtLower Var" << endl;
                        triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], 1));
                        //A_basis(basis_i, var_index_map[var.getId()]) = -1;
                        ++basis_i;
                    }
                }

                //cin.ignore();

                A_basis.setFromTriplets(triplets.begin(), triplets.end());

                vars_basis.end();
                constrs_basis.end();
            }

            {
                IloNumArray cplex_optLP(cplex.getEnv());
                cplex.getValues(vars, cplex_optLP);

                for (size_t i = 0; i < n; ++i) {
                    if (cplex_optLP[i] != 0)
                        optLP.insertBack(i) = cplex_optLP[i];
                }

                cplex_optLP.end();
            }

            A_basis = adjust_normal_cone(A_basis.transpose(), c);

            return std::move(std::make_tuple(A_basis, c, optLP));
        }
    };
}
