
#include <iostream>
#include <random>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "conic_cuts_Eigen.h"
#include "mpir_eigen_supp.h"

namespace ILPCuts {
    template<typename intT = long long int, typename ratT = double>
    class CPlexAdapter :public ConicCutsGenerator<intT, ratT> {
    public:

        using base = ConicCutsGenerator<intT, ratT>;
        /*using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<intT, ratT>;

        using int_type = typename CutsGen::int_type;
        using rat_type = typename CutsGen::rational_type;
        using rat_vec = typename CutsGen::rat_vec;
        using int_vec = typename CutsGen::int_vec;
        using rat_mat = typename CutsGen::rat_mat;
        using int_mat = typename CutsGen::int_mat;*/

        using typename base::int_type;
        using typename base::rational_type;
        using typename base::is_mat;
        using typename base::rs_mat;
        using typename base::is_vec;
        using typename base::rs_vec;

        CPlexAdapter(std::ostream* outs, std::ostream* errs)
            :base(outs, errs) {}

        CPlexAdapter() :base() {}

        void add_cuts(IloCplex cplex, IloNumVarArray vars, const std::vector<is_vec> & cuts_Eigen) {
            IloEnv env = cplex.getEnv();

            IloConstraintArray ilo_cuts(env);

            size_t n = vars.getSize();

            //std::map<IloInt, size_t> var_index_map;
            //for (size_t i = 0; i < n; ++i) {
            //    var_index_map[vars[i].getId()] = i;
            //}
            /*for (size_t i = 0; i < cuts_Eigen.size(); ++i) {
                base::out() << cuts_Eigen[i] << endl;
            }*/

            for (size_t i = 0; i < cuts_Eigen.size(); ++i) {
                IloExpr lhs = IloExpr(env);
                IloNum rhs = 0;

                for (typename is_vec::InnerIterator it(cuts_Eigen[i]); it; ++it) {
                    if (it.index() < n)
                        lhs += vars[it.index()] * it.value().get_si();
                    else
                        rhs = it.value().get_si();
                }

                /*for (size_t j = 0; j+1 < n; ++j) {
                    lhs += vars[j] * cuts_Eigen[i][j].get_si();
                }*/
                //IloNum rhs = cuts_Eigen[i].coeff(n).get_si();
                ilo_cuts.add( IloRange(env, lhs, rhs) );
            }

            cplex.addUserCuts(ilo_cuts);

            ilo_cuts.end();
        }

        void add_cuts(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars) {
            is_mat A;
            rs_vec optLP;

            tie(A, optLP) = transform_to_corner_polyhedron(cplex, model, obj, constrs, vars);

            std::mt19937 engine;
            typename base::template ConvUniform<std::mt19937> ds(engine, A);

            //typename base::SomeDirection ds(A);

            std::vector<is_vec> cuts_Eigen = base::generate_cuts(A, optLP, ds, 0);
            add_cuts(cplex, vars, cuts_Eigen);
        }

        std::tuple<is_mat, rs_vec>
            transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars) {
            size_t n = vars.getSize();

            typename base::is_mat A_basis(n,n);
            rs_vec optLP(n), c(n);

            std::map<IloInt, size_t> var_index_map;
            for (size_t i = 0; i < n; ++i) {
                var_index_map[vars[i].getId()] = i;
            }

            {
                IloExpr::LinearIterator iter = obj.getLinearIterator();
                while (iter.ok()) {
                    if (obj.getSense() == IloObjective::Sense::Maximize)
                        c.insertBack(var_index_map[iter.getVar().getId()]) = iter.getCoef();
                    else
                        c.insertBack(var_index_map[iter.getVar().getId()]) = -iter.getCoef();
                    ++iter;
                }
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
                    if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], iter.getCoef()));
                            //A_basis(basis_i, var_index_map[iter.getVar().getId()]) = iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                    else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            triplets.push_back(Trilpl(basis_i, var_index_map[iter.getVar().getId()], -iter.getCoef()));
                            //A_basis(basis_i, var_index_map[iter.getVar().getId()]) = -iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                }

                for (size_t i = 0; i < vars.getSize(); ++i) {
                    IloNumVar var = vars[i];
                    if (vars_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], 1));
                        //A_basis(basis_i, var_index_map[var.getId()]) = 1;
                        ++basis_i;
                    }
                    else if (vars_basis[i] == IloCplex::BasisStatus::AtLower) {
                        triplets.push_back(Trilpl(basis_i, var_index_map[var.getId()], -1));
                        //A_basis(basis_i, var_index_map[var.getId()]) = -1;
                        ++basis_i;
                    }
                }

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
                    //optLP[i] = cplex_optLP[i];

                cplex_optLP.end();
            }

            return std::move(std::make_tuple(A_basis.transpose(), optLP));
        }

        /*std::tuple<is_mat, rs_vec>
            transform_to_corner_polyhedron(IloCplex cplex, IloModel model, IloObjective obj, IloRangeArray constrs, IloNumVarArray vars) {
            size_t n = vars.getSize();

            typename base::int_mat A_basis = int_mat::Zero(n, n);
            rat_vec optLP = rat_vec::Zero(n);
            rat_vec c = rat_vec::Zero(n);

            std::map<IloInt, size_t> var_index_map;
            for (size_t i = 0; i < n; ++i) {
                var_index_map[vars[i].getId()] = i;
            }

            {
                IloExpr::LinearIterator iter = obj.getLinearIterator();
                while (iter.ok()) {
                    if (obj.getSense() == IloObjective::Sense::Maximize)
                        c[var_index_map[iter.getVar().getId()]] = iter.getCoef();
                    else
                        c[var_index_map[iter.getVar().getId()]] = -iter.getCoef();
                    ++iter;
                }
            }

            {
                IloCplex::BasisStatusArray vars_basis(cplex.getEnv());
                IloCplex::BasisStatusArray constrs_basis(cplex.getEnv());
                cplex.getBasisStatuses(vars_basis, vars, constrs_basis, constrs);
                size_t basis_i = 0;

                for (size_t i = 0; i < constrs.getSize(); ++i) {
                    IloRange constr = constrs[i];
                    if (constrs_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            A_basis(basis_i, var_index_map[iter.getVar().getId()]) = iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                    else if (constrs_basis[i] == IloCplex::BasisStatus::AtLower) {
                        IloExpr::LinearIterator iter = constr.getLinearIterator();
                        while (iter.ok()) {
                            A_basis(basis_i, var_index_map[iter.getVar().getId()]) = -iter.getCoef();
                            ++iter;
                        }
                        ++basis_i;
                    }
                }

                for (size_t i = 0; i < vars.getSize(); ++i) {
                    IloNumVar var = vars[i];
                    if (vars_basis[i] == IloCplex::BasisStatus::AtUpper) {
                        A_basis(basis_i, var_index_map[var.getId()]) = 1;
                        ++basis_i;
                    }
                    else if (vars_basis[i] == IloCplex::BasisStatus::AtLower) {
                        A_basis(basis_i, var_index_map[var.getId()]) = -1;
                        ++basis_i;
                    }
                }

                vars_basis.end();
                constrs_basis.end();
            }

            {
                IloNumArray cplex_optLP(cplex.getEnv());
                cplex.getValues(vars, cplex_optLP);

                for (size_t i = 0; i < n; ++i)
                    optLP[i] = cplex_optLP[i];

                cplex_optLP.end();
            }

            return std::move(std::make_tuple(A_basis.transpose(), optLP));
        }*/
    };
}

//ILOUSERCUTCALLBACK4(myCutsCallback, IloCplex, cplex, IloObjective, obj, IloNumVarArray, vars, IloRangeArray, rangs) {
//    IloEnv env = getEnv();
//    IloModel mod = getModel();
//
//    if (getNnodes() == 0) {
//        //env.out() << "--------------My cut Callback----------------" << std::endl;
//
//        using CutsGen = ILPCuts::CPlexAdapter<mpz_class, mpz_class>;
//        CutsGen cuts_gen;
//
//        typename CutsGen::int_mat A;
//        typename CutsGen::rat_vec optLP;
//
//        tie(A, optLP) = cuts_gen.transform_to_corner_polyhedron(cplex, mod, obj, rangs, vars);
//
//        //env.out() << "---------------------------------------------" << std::endl;
//    }
//
//    //{A, optLP} = cuts_gen.transform_to_corner_polyhedron(cplex, model);
//
//    /*IloIterator<IloExpr> citer = IloIterator<IloExpr>(env);
//    IloModel::Iterator iter = IloModel::Iterator(mod);
//    IloConstraintArray carr(env);
//    IloNumVarArray varr(env);
//    IloExprArray earr(env);
//    while (iter.ok()) {
//        if ((*iter).isConstraint())
//            carr.add((*iter).asConstraint());
//        else if ((*iter).isVariable())
//            varr.add((*iter).asVariable());
//        else if ((*iter).isNumExpr())
//            earr.add((*iter).asNumExpr());
//        else if ((*iter).isIntExpr())
//            earr.add((*iter).asIntExpr());
//        ++iter;
//    }*/
//
//    /*env.out() << "=====variables==========" << endl;
//    for (size_t i = 0; i < varr.getSize(); ++i) {
//        env.out() << varr[i].getName() << ", ";
//    }
//    env.out() << endl;
//
//    env.out() << "=====constraints==========" << endl;
//    for (size_t i = 0; i < carr.getSize(); ++i) {
//        env.out() << carr[i].getName() << ", ";
//    }
//    env.out() << endl;*/
//
//    //env.out() << "=====constraints 2.0==========" << endl;
//    //size_t expr_count = 0;
//    //while (citer.ok()) {
//    //    if ((*citer).isNumExpr())
//    //    {
//    //        //env.out() << (*citer).asNumExpr() << ", ";
//    //        ++expr_count;
//    //    }
//    //    else if ((*citer).isIntExpr()) {
//    //        //env.out() << (*citer).asIntExpr() << ", ";
//    //        ++expr_count;
//    //    }
//    //        
//    //    ++citer;
//    //    
//    //}
//
//    //env.out()<< expr_count << '\n';
//    /*
//    env.out() << "=====exprs==========" << endl;
//    for (size_t i = 0; i < earr.getSize(); ++i) {
//        env.out() << earr[i].getName() << ", ";
//    }
//    
//    */
//
//    /*carr.end();
//    varr.end();
//    earr.end();*/
//    //mod.end();
//}

ILOSIMPLEXCALLBACK0(SimplexCollback) {
    IloEnv env = getEnv();

}

int main(int argc, char** argv)
{
    IloEnv env;
    try {
        IloModel master_model(env);
        IloModel relaxed_model(env);
        IloCplex cplex(env);
        IloObjective   obj;
        IloNumVarArray vars(env);
        IloRangeArray  rngs(env);

        const char* datadir = (argc >= 2) ? argv[1] : "../problems/cplex";
        stringstream filename;
        filename << datadir << "/" << "noswot.mps";
        env.out() << "reading " << filename.str() << endl;
        cplex.importModel(master_model, filename.str().c_str(), obj, vars, rngs);

        relaxed_model.add(master_model);
        for (size_t i = 0; i < vars.getSize(); ++i) {
            relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
        }

        /*cplex.extract(master_model);
        env.out() << "is MIP: " << cplex.isMIP() << endl;
        cplex.extract(relaxed_model);
        env.out() << "is MIP: " << cplex.isMIP() << endl;*/

        /*env.out() << "extracting master model..." << endl;
        cplex.extract(master_model);*/

        cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
        cplex.extract(relaxed_model);
        env.out() << "solving relaxed model ...\n";

        //cplex.use(myCutsCallback(env, cplex, obj, vars, rngs));

        if (!cplex.solve()) {
            env.error() << "Failed to optimize." << endl;
            throw(-1);
        }

        env.out() << "solution status is " << cplex.getStatus() << endl;
        env.out() << "solution value  is " << cplex.getObjValue() << endl;

        using CutsGen = ILPCuts::CPlexAdapter<mpz_class, mpq_class>;
        CutsGen cuts_gen(&env.out(), &env.out());

        env.out() << "cutting generation:" << endl;
        cuts_gen.add_cuts(cplex, relaxed_model, obj, rngs, vars);

        cin.ignore();

        cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
        cplex.extract(master_model);
        env.out() << "solving relaxed model ...\n";

        if (!cplex.solve()) {
            env.error() << "Failed to optimize." << endl;
            throw(-1);
        }

        env.out() << "solution status is " << cplex.getStatus() << endl;
        env.out() << "solution value  is " << cplex.getObjValue() << endl;
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }

    env.end();

    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
