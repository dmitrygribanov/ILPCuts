
#include <iostream>
#include <random>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "cplex_adapter.h"
#include "conic_cuts_Eigen.h"
#include "conic_cuts_Eigen_optimized.h"
#include "../headers/argparse/argparse.hpp"
using namespace ILPCuts;

using Strategies = typename CutsGenStrategies<mpz_class, mpf_class>;
using CutsGenAdapter = typename CPlexAdapter<mpz_class, mpf_class>;
using CutsGen_Adj = typename ConicCutsGenerator<mpz_class, mpf_class>;
using CutsGen_LU = typename ConicCutsGeneratorOpt<mpz_class, mpf_class, mpz_class, mpf_class>;

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

//enum class SolverTypeEnum {Adj, LU};

void run_cplex_experiment(
	std::string input_path,
	double tolerance,
	double err_tolerance,
	size_t cone_limit_per_direction,
	std::string output_path,
	size_t rand_directions,
	bool some_direction,
	size_t threads,
	bool baseline,
	bool lazy_constraints_on,
	size_t iterations_number,
	std::string solver_type) {
	//TO DO: дописать сохраниение в формат csv
	//create CPLEX envirement
	IloEnv env;

	//set global presision for long arithmetics
	mpf_set_default_prec(200);

	//create main model object
	IloModel master_model(env);

	//create relaxed model object
	IloModel relaxed_model(env);

	//create solver object
	IloCplex cplex_main(env);

	//solver object to solve relaxed model
	IloCplex cplex_relaxed(env);

	//create objective function object
	IloObjective   obj;

	//object for variables
	IloNumVarArray vars(env);

	//object for inequalities
	IloRangeArray  rngs(env);

	env.out() << "reading " << input_path << endl;

	//import main problem from file
	cplex_main.importModel(master_model, input_path.c_str(), obj, vars, rngs);

	//import data to the relaxed model
	relaxed_model.add(master_model);

	//make all variables of the relaxed model to be rational
	for (size_t i = 0; i < vars.getSize(); ++i) {
		relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
	}

	//setting parameters to solve relaxed model
	cplex_relaxed.setParam(IloCplex::Param::MIP::Interval, 1000);
	cplex_relaxed.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	//put relaxed model to the solver
	cplex_relaxed.extract(relaxed_model);
	env.out() << "solving relaxed model ...\n";

	//run solver
	if (!cplex_relaxed.solve()) {
		env.error() << "Failed to optimize." << endl;
		throw(-1);
	}
	env.out() << "solution status is " << cplex_relaxed.getStatus() << endl;
	env.out() << "solution value  is " << cplex_relaxed.getObjValue() << endl;

	//return to the original model
	cplex_main.extract(master_model);

	//createing cuts generator object
	CutsGenAdapter cuts_gen_adapter;
	if (solver_type == "adj")
		cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_Adj>(&env.out(), &env.error());
	else if (solver_type == "lu")
		cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_LU>(&env.out(), &env.error());
	else {
		env.error() << "Wrong value of the CLI parameter -cg(--cgen): " << solver_type << std::endl;
		throw(-1);
	}
	auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

	//initialization of cuts generator and setting og its parameters
	std::mt19937 engine;
	cuts_gen_ptr->Init();
	cuts_gen_ptr->tolerance = tolerance;
	cuts_gen_ptr->err_tolerance = err_tolerance;
	cuts_gen_ptr->cone_limit_per_direction = cone_limit_per_direction;
	cuts_gen_ptr->threads = threads;
	cuts_gen_adapter.lazy_constraints_on = lazy_constraints_on;
	cuts_gen_ptr->rgen = &engine;

	//set strategies for cutting generation
	if (rand_directions > 0) {
		using RandStrategy = typename Strategies::ConvUniform<std::mt19937>;
		cuts_gen_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
	}
	if (some_direction)
		cuts_gen_ptr->direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());

	if (baseline == false) {
		env.out() << "cutting generation started:" << endl;

		for (size_t iteration_idx = 0; iteration_idx < iterations_number; ++iteration_idx) {
			env.out() << "root node boosting iteration " << iteration_idx << " is started :" << endl;
			//cin.ignore();

			//generate cuts from relaxed model
			env.out() << "constraints generation ..." << endl;
			cuts_gen_adapter.generate_cuts(cplex_relaxed, relaxed_model, obj, rngs, vars);

			env.out() << "adding constraints to relaxed model ..." << endl;
			cuts_gen_adapter.add_cuts_as_constraints(cplex_relaxed, relaxed_model, rngs, vars);

			env.out() << "adding cuts to original model ..." << endl;
			cuts_gen_adapter.add_cuts(cplex_main, vars);

			//put relaxed model to the solver
			cplex_relaxed.extract(relaxed_model);
			env.out() << "solving boosted model ...\n";

			//run solver
			if (!cplex_relaxed.solve()) {
				env.error() << "Failed to optimize." << endl;
				throw(-1);
			}

			env.out() << "solution status is " << cplex_relaxed.getStatus() << endl;
			env.out() << "solution value  is " << cplex_relaxed.getObjValue() << endl;
		}
		if (iterations_number > 0)
			env.out() << "root node boosting finished." << endl;

		cuts_gen_adapter.generate_cuts(cplex_relaxed, relaxed_model, obj, rngs, vars);
		cuts_gen_adapter.add_cuts(cplex_main, vars);

		env.out() << "cutting generation finished:" << endl;
	}

	//setting parameters for solving in integer variables
	cplex_main.setParam(IloCplex::Param::Threads, threads);
	/*cplex_main.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
	cplex_main.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);*/

	//cplex_main.setParam(IloCplex::Param::MIP::Interval, 1000);
	//cplex_main.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	//solve the ILP model
	if (!cplex_main.solve()) {
		env.error() << "Failed to optimize." << endl;
		throw(-1);
	}

	env.out() << "solution status is " << cplex_main.getStatus() << endl;
	env.out() << "solution value  is " << cplex_main.getObjValue() << endl;

	env.end();
}

int main(int argc, char** argv)
{
	argparse::ArgumentParser program("IPLCutsCLI");

	program.add_argument("-i", "--input")
		.default_value(std::string("-"))
		.required()
		.help("specify the input file.");

	program.add_argument("-o", "--output")
		.default_value(std::string(""))
		.help("specify the output file. (not implemented yet)");

	program.add_argument("-cg", "--cgen")
		.default_value(std::string("adj"))
		.help("specify the output file. (not implemented yet)");

	program.add_argument("-t", "--tolerance")
		.default_value(0.0000001)
		.help("ILPCuts algorithm tolerance")
		.scan<'g', double>();

	program.add_argument("-e", "--err_tolerance")
		.default_value(0.0001)
		.help("ILPCuts algorithm error tolerance")
		.scan<'g', double>();

	program.add_argument("-c", "--cone_limit")
		.default_value(size_t(1))
		.help("cone limit per direction")
		.scan<'i', size_t>();

	program.add_argument("-r", "--rand_strategy")
		.default_value(size_t(0))
		.help("cone limit in random cone strategy")
		.scan<'i', size_t>();

	program.add_argument("-s", "--some_strategy")
		.default_value(false)
		.help("use of some cone strategy ")
		.implicit_value(true);

	program.add_argument("-w", "--workers", "--threads")
		.default_value(size_t(1))
		.help("number of threads")
		.scan<'i', size_t>();

	program.add_argument("-b", "--baseline")
		.default_value(false)
		.help("use baseline solver only")
		.implicit_value(true);

	program.add_argument("-lz", "--lazy")
		.default_value(false)
		.help("turns on addition of CPlex lazy constraints")
		.implicit_value(true);

	program.add_argument("-it", "--iterations")
		.default_value(size_t(0))
		.help("number of iterations in the root-node")
		.scan<'i', size_t>();

	try {
		program.parse_args(argc, argv);

		run_cplex_experiment(
			program.get<std::string>("--input"),
			program.get<double>("--tolerance"),
			program.get<double>("--err_tolerance"),
			program.get<size_t>("--cone_limit"),
			program.get<std::string>("--output"),
			program.get<size_t>("--rand_strategy"),
			program.get<bool>("--some_strategy"),
			program.get<size_t>("--threads"),
			program.get<bool>("--baseline"),
			program.get<bool>("--lazy"),
			program.get<size_t>("--iterations"),
			program.get<std::string>("--cgen")
		);
	}
	catch (const std::exception& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << program;
		std::exit(1);
	}

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
