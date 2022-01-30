
#include <iostream>
#include <random>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "cplex_adapter.h"
#include "conic_cuts_Eigen.h"
#include "mpir_eigen_supp.h"

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

using namespace ILPCuts;

int main(int argc, char** argv)
{
    //create CPLEX envirement
    IloEnv env;

    //set global presision for long arithmetics
    mpf_set_default_prec(200);


    try {
        //create main model object
        IloModel master_model(env);

        //create relaxed model object
        IloModel relaxed_model(env);

        //create solver object
        IloCplex cplex(env);

        //create objective function object
        IloObjective   obj;

        //object for variables
        IloNumVarArray vars(env);

        //object for inequalities
        IloRangeArray  rngs(env);

        //defoult datav=bse dir
        const char* datadir = (argc >= 2) ? argv[1] : "../problems/cplex";
        
        stringstream filename;
        filename << datadir << "/" << "sorrell3.mps.gz";
        //filename << datadir << "/" << "noswot.mps";
        //filename << datadir << "/" << "bnatt500.mps.gz";
        //filename << datadir << "/" << "splice1k1.mps.gz";
        //filename << datadir << "/" << "fhnw-binpack4-48.mps.gz";
        //filename << datadir << "/" << "hanoi5.mps.gz";
        //filename << datadir << "/" << "b1c1s1.mps.gz";
        
        
        env.out() << "reading " << filename.str() << endl;

        //import main problem from file
        cplex.importModel(master_model, filename.str().c_str(), obj, vars, rngs);

        //import data to the relaxed model
        relaxed_model.add(master_model);

        //make all variables of the relaxed model to be rational
        for (size_t i = 0; i < vars.getSize(); ++i) {
            relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
        }

        //setting parameters to solve relaxed model
        cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
        cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

        //put relaxed model to the solver
        cplex.extract(relaxed_model);
        env.out() << "solving relaxed model ...\n";

        //cplex.use(myCutsCallback(env, cplex, obj, vars, rngs));

        //run solver
        if (!cplex.solve()) {
            env.error() << "Failed to optimize." << endl;
            throw(-1);
        }
        env.out() << "solution status is " << cplex.getStatus() << endl;
        env.out() << "solution value  is " << cplex.getObjValue() << endl;

        //createing cuts generator object
        using CutsGen = CPlexAdapter<ConicCutsGenerator<mpz_class, mpf_class> >;
        CutsGen cuts_gen(&env.out(), &env.out());

        //initialization of cuts generator and setting og its parameters
        cuts_gen.Init();
        cuts_gen.tolerance = 0.000001;
        cuts_gen.err_tolerance = 0.001;
        cuts_gen.cone_limit_per_direction = 6;

        //set strategies for cutting generation
        std::mt19937 engine;
        using RandStrategy = typename CutsGen::template ConvUniform<std::mt19937>;
        //cuts_gen.direction_strategies.push_back(std::make_unique<RandStrategy>(engine, 3));
        //cuts_gen.direction_strategies.push_back(std::make_unique<typename CutsGen::SomeDirection>());

        //generate cuts from relaxed model
        env.out() << "cutting generation:" << endl;
        cuts_gen.generate_cuts(cplex, relaxed_model, obj, rngs, vars);

        //put the original model to the solver
        cplex.extract(master_model);

        //adding our cuts to the solver
        cuts_gen.add_cuts(cplex, vars);


        //setting parameters for solving in integer variables
        cplex.setParam(IloCplex::Param::Threads, 7);
        /*cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
        cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);*/
        //cin.ignore();

        //cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
        //cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

        //solve the ILP model
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
    catch (std::exception & ex) {
        cerr << ex.what() << endl;
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
