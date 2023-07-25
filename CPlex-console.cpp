#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

// #define ILPCUTS_PARALLEL_OFF

#include "cone_coverer_adj.h"
#include "cone_coverer_LU.h"
#include "cone_coverer_Zp.h"
#include "cone_coverer_Wd.h"
#include "cplex_adapter.h"
#include "cuts_gen_submat.h"
#include "dependencies/argparse/argparse.hpp"
#include "ILPCuts_console_utils.h"

using namespace ILPCuts;
using namespace ModuleCPlex;

void check_input_params(
	logger &log,
	double tolerance,
	double err_tolerance,
	size_t cone_limit_per_direction,
	size_t threads_global,
	string is_lazy_param,
	int node_lim,
	std::string solver_type,
	size_t cprime,
	string cplex_cuts_param,
	string scaling,
	string model_cuts_switch,
	string cuts_fun_type)
{
	if (tolerance <= 0.0)
	{
		log.out() << "Wrong value of the CLI parameter -t(--tolerance): "
				  << tolerance << std::endl;
		exit(0);
	}

	if (err_tolerance <= 0.0)
	{
		log.out() << "Wrong value of the CLI parameter -e(--err_tolerance): "
				  << err_tolerance << std::endl;
		exit(0);
	}

	if (cone_limit_per_direction == 0)
	{
		log.out() << "Wrong value of the CLI parameter -c(--cone_limit): "
				  << cone_limit_per_direction << std::endl;
		exit(0);
	}

	if (threads_global == 0)
	{
		log.out() << "Wrong value of the CLI parameter -w(--workers, --threads): "
				  << threads_global << std::endl;
		exit(0);
	}

	if (cprime < 2)
	{
		log.out() << "Wrong value of the CLI parameter --cprime: "
				  << cprime << std::endl;
		exit(0);
	}

	if (solver_type != "adj" && solver_type != "lu" && solver_type != "zp" && solver_type != "wd")
	{
		log.out() << "Wrong value of the CLI parameter -cg(--cgen): "
				  << solver_type << std::endl;
		exit(0);
	}

	if (is_lazy_param != "on" && is_lazy_param != "off" && is_lazy_param != "hybrid")
	{
		log.out() << "Wrong value of the CLI parameter --lazy: "
				  << is_lazy_param << endl;
		exit(0);
	}

	// if (iterations_number == 0)
	// {
	// 	log.out() << "Wrong value of the CLI parameter -it(--iterations): "
	// 			  << iterations_number << std::endl;
	// 	exit(0);
	// }

	if (model_cuts_switch != "on" && model_cuts_switch != "off")
	{
		log.out() << "Wrong value of the CLI parameter --model_cuts: "
				  << model_cuts_switch << std::endl;
		exit(0);
	}

	if (node_lim < -1)
	{
		log.out() << "Wrong value of the CLI parameter --node_lim: "
				  << node_lim << std::endl;
		exit(0);
	}

	if (scaling != "on" && scaling != "off")
	{
		log.out() << "Wrong value of the CLI parameter --scaling: "
				  << scaling << std::endl;
		exit(0);
	}

	if (cplex_cut_types.check_name(cplex_cuts_param))
	{
		log.out() << "Wrong value of the CLI parameter --cplex_cuts_param: "
				  << cplex_cuts_param << std::endl;
		exit(0);
	}

	if (cuts_fun_type != "split" && cuts_fun_type != "mir" && cuts_fun_type != "hermite" && cuts_fun_type != "pureint")
	{
		log.out() << "Wrong value of the CLI parameter --cuts_fun: "
				  << cuts_fun_type << std::endl;
		exit(0);
	}
}

void run_cplex_experiment(
	logger &log,
	string input_path,
	string cuts_output_path,
	double tolerance,
	double err_tolerance,
	size_t cone_limit_per_direction,
	size_t rand_directions,
	bool some_direction,
	size_t threads_cuts_gen,
	size_t threads_cplex,
	string is_lazy_param,
	size_t iterations_number,
	int node_lim,
	string solver_type,
	size_t cprime,
	string cplex_cuts_param,
	string scaling,
	string model_cuts_switch,
	size_t min_dim,
	size_t max_dim,
	double cplex_tl,
	double cuts_tl,
	size_t cuts_generation_submat_size_limit,
	string cuts_fun_type,
	size_t row_trials_num)
{
	// create CPLEX envirement
	IloEnv env;
	env.setOut(log.out());
	env.setError(log.out());

	// set global presision for long arithmetics
	mpf_set_default_prec(256);

	// create main model object
	IloModel master_model(env);

	// create relaxed model object
	IloModel relaxed_model(env);

	// create solver object for the main MILP
	IloCplex cplex_master(env);

	// solver object to solve the relaxed model
	IloCplex cplex_relaxed(env);

	// create objective function object
	IloObjective obj;

	// object for variables
	IloNumVarArray vars(env);

	// object for inequalities
	IloRangeArray rngs(env);

	env.out() << "reading the problem..." << input_path << endl;

	// import main problem from file
	cplex_master.importModel(master_model, input_path.c_str(), obj, vars, rngs);

	// env.out() << print_IloRangeArray(rngs) << endl;

	// checking dimensional requrements for the problem
	if (vars.getSize() < min_dim)
	{
		log.out() << "\nDimension of the problem (" << vars.getSize()
				  << ") is less than min_dim (" << min_dim << ") parameter,\n"
				  << "terminate..." << endl;
		exit(0);
	}
	if (max_dim > 0 && vars.getSize() > max_dim)
	{
		log.out() << "\nDimension of the problem (" << vars.getSize()
				  << ") is greater than max_dim (" << max_dim << ") parameter,\n"
				  << "terminate..." << endl;
		exit(0);
	}

	// import data to the relaxed model
	relaxed_model.add(master_model);

	// make all variables of the relaxed model to be rational
	for (size_t i = 0; i < vars.getSize(); ++i)
	{
		relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
	}

	// setting parameters to solve relaxed model
	// cplex_relaxed.setParam(IloCplex::Param::MIP::Interval, 1000);
	cplex_relaxed.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	// createing adapter object
	CPlexAdapter cuts_gen_adapter;

	if (cuts_fun_type == "split")
	{
		cuts_gen_adapter.generator_ptr = std::make_unique<SplitCutsGen_RWDS>(log);
	}
	else if (cuts_fun_type == "hermite")
	{
		cuts_gen_adapter.generator_ptr = std::make_unique<HermiteCutsGen_RWDS>(log);
	}
	else if (cuts_fun_type == "pureint")
	{
		cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_CC>(log);
	}
	else if (cuts_fun_type == "mir")
	{
		cuts_gen_adapter.generator_ptr = std::make_unique<MIRCutsGen_RWDS>(log);
	}
	else
	{
		log.out() << "Wrong value of the CLI parameter --cuts_fun: "
				  << cuts_fun_type << std::endl;
		exit(0);
	}
	auto cuts_gen_ptr = static_cast<CutsGen_CC *>(cuts_gen_adapter.generator_ptr.get());

	// creating and wrapping of a cuts generator object
	if (solver_type == "adj")
		cuts_gen_ptr->cone_coverer_ptr = std::make_unique<ConeCoverer_adj>(log);
	else if (solver_type == "lu")
		cuts_gen_ptr->cone_coverer_ptr = std::make_unique<ConeCoverer_LU>(log);
	else if (solver_type == "zp")
		cuts_gen_ptr->cone_coverer_ptr = std::make_unique<ConeCoverer_Zp>(log);
	else if (solver_type == "wd")
		cuts_gen_ptr->cone_coverer_ptr = std::make_unique<ConeCoverer_Wd>(log);

	auto cone_coverer_ptr = cuts_gen_ptr->cone_coverer_ptr.get();

	// initialization of cuts generator and setting parameters
	std::mt19937 engine;
	cuts_gen_ptr->set_params(tolerance, err_tolerance, threads_cuts_gen, &engine);
	if (cuts_fun_type == "split")
	{
		static_cast<SplitCutsGen_RWDS *>(cuts_gen_ptr)->submat_size_fun.maxSize = cuts_generation_submat_size_limit;
		static_cast<SplitCutsGen_RWDS *>(cuts_gen_ptr)->submat_gen.trials_num = row_trials_num;
	}
	else if (cuts_fun_type == "hermite")
	{
		static_cast<HermiteCutsGen_RWDS *>(cuts_gen_ptr)->submat_size_fun.maxSize = cuts_generation_submat_size_limit;
		static_cast<HermiteCutsGen_RWDS *>(cuts_gen_ptr)->submat_gen.trials_num = row_trials_num;
	}
	else if (cuts_fun_type == "mir")
	{
		static_cast<MIRCutsGen_RWDS *>(cuts_gen_ptr)->submat_size_fun.maxSize = cuts_generation_submat_size_limit;
		static_cast<MIRCutsGen_RWDS *>(cuts_gen_ptr)->submat_gen.trials_num = row_trials_num;
	}

	cone_coverer_ptr->cone_limit_per_direction = cone_limit_per_direction;
	cuts_gen_adapter.set_lazy(is_lazy_param);
	if (solver_type == "zp")
	{
		static_cast<ConeCoverer_Zp *>(cone_coverer_ptr)->max_conservative_prime = cprime;
	}
	else if (solver_type == "wd")
	{
		static_cast<ConeCoverer_Wd *>(cone_coverer_ptr)->max_conservative_prime = cprime;
	}

	// set strategies for cutting generation
	if (rand_directions > 0)
	{
		using RandStrategy = typename DirectionStrategies::ConvUniform<std::mt19937>;
		cone_coverer_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
	}
	if (some_direction)
		cone_coverer_ptr->direction_strategies.push_back(std::make_unique<DirectionStrategies::SomeDirection>());

	// setting parameters for solving in integer variables
	cplex_master.setParam(IloCplex::Param::Threads, threads_cplex);

	if (cplex_tl > 0.0)
	{
		cplex_master.setParam(IloCplex::Param::TimeLimit, cplex_tl);
		// cplex_master.setParam(IloCplex::Param::DetTimeLimit, 0);
	}

	if (node_lim != -1)
		cplex_master.setParam(IloCplex::Param::MIP::Limits::Nodes, node_lim);

	if (scaling == "off")
	{
		cplex_master.setParam(IloCplex::Param::Read::Scale, -1);
		// cplex_master.setParam(IloCplex::Param::MIP::SubMIP::Scale, -1);
	}

	// setting cuts type for cplex_master
	cplex_cut_types.set_param(cplex_master, cplex_cuts_param);

	// cplex_master.setParam(IloCplex::Param::MIP::Interval, 1000);
	// cplex_master.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	double originalRelaxedObjVal = root_node_boosting_console_subprogram(cplex_master, cplex_relaxed,
																		 obj, rngs, vars,
																		 cuts_gen_adapter,
																		 cuts_tl, iterations_number, model_cuts_switch, cuts_output_path);

	double boostedRelaxedObjVal = cplex_relaxed.getObjValue();

	// solve the ILP model
	if (cplex_master.solve())
	{
		env.out() << "ILP solution value: " << cplex_master.getObjValue() << "\n"
				  << "ILP cutoff: " << cplex_master.getCutoff() << endl;

		if (iterations_number != 0)
		{
			env.out()
				<< "Duality gap decrease (ILP): "
				<< (originalRelaxedObjVal - boostedRelaxedObjVal) / (originalRelaxedObjVal - cplex_master.getObjValue()) << "\n"
				<< "Duality gap decrease (cutoff): "
				<< (originalRelaxedObjVal - boostedRelaxedObjVal) / (originalRelaxedObjVal - cplex_master.getCutoff()) << "\n"
				<< endl;
		}
	}
	else
	{
		env.out() << "Integer optimal solution is not found...\n"
				  << "Solution status is " << cplex_master.getStatus() << "\n"
				  << endl;
	}

	env.out()
		<< "LP best solution value: " << cplex_master.getBestObjValue() << endl;

	if (iterations_number != 0)
	{
		env.out() << "Duality gap decrease (LP): "
				  << (originalRelaxedObjVal - boostedRelaxedObjVal) / (originalRelaxedObjVal - cplex_master.getBestObjValue()) << endl;
	}

	env.end();
}

int main(int argc, char **argv)
{
	argparse::ArgumentParser program("IPLCutsCLI");

	program.add_argument("-i", "--input")
		.default_value(std::string())
		.required()
		.help("specify the input file.");

	program.add_argument("-o", "--output")
		.default_value(std::string("log.txt"))
		.help("output log file; default is \"log.txt\"");

	program.add_argument("--cuts_output")
		.default_value(std::string("cuts_export.txt"))
		.help("cuts export file; default is \"cuts_export.txt\"");

	program.add_argument("--cgen")
		.default_value(std::string("zp"))
		.help("type of cutting generation algorithm; can be \"adj\", \"lu\", \"zp\" or \"wd\"");

	program.add_argument("-t", "--tolerance")
		.default_value(0.000001)
		.help("ILPCuts algorithm tolerance")
		.scan<'g', double>();

	program.add_argument("--cprime")
		.default_value(size_t(100000))
		.help("max conservative prime for ILPCuts with Zp")
		.scan<'i', size_t>();

	program.add_argument("-e", "--err_tolerance")
		.default_value(0.0001)
		.help("ILPCuts algorithm error tolerance")
		.scan<'g', double>();

	program.add_argument("-c", "--cone_limit")
		.default_value(size_t(1))
		.help("cone limit per direction")
		.scan<'i', size_t>();

	program.add_argument("--submat_limit")
		.default_value(size_t(1024))
		.help("square submatrix size limit")
		.scan<'i', size_t>();

	program.add_argument("-r", "--rand_strategy")
		.default_value(size_t(0))
		.help("number of directions in random cone strategy")
		.scan<'i', size_t>();

	program.add_argument("--sum_strat")
		.default_value(false)
		.help("add \"sum\" cone strategy")
		.implicit_value(true);

	program.add_argument("-w", "--workers", "--threads")
		.default_value(size_t(1))
		.help("global number of threads")
		.scan<'i', size_t>();

	program.add_argument("--wcgn", "--workers_cuts", "--threads_cuts")
		.default_value(size_t(0))
		.help("number of threads for cutting generation algorithm")
		.scan<'i', size_t>();

	program.add_argument("--wcpl", "--workers_cplex", "--threads_cplex")
		.default_value(size_t(0))
		.help("number of threads for CPlex solver")
		.scan<'i', size_t>();

	program.add_argument("-b", "--baseline")
		.default_value(false)
		.help("use baseline solver only")
		.implicit_value(true);

	program.add_argument("--lazy")
		.default_value(string("off"))
		.help("manages CPlex lazy constraints; can be \"on\", \"off\" or \"hybrid\"");

	program.add_argument("--iters")
		.default_value(size_t(1))
		.help("number of boosting iterations in the root-node")
		.scan<'i', size_t>();

	program.add_argument("--node_lim")
		.default_value(int(-1))
		.help("maximal number of nodes in CPlex B&C-algorithm (-1 = default CPlex behavor)")
		.scan<'i', int>();

	program.add_argument("--cplex_cuts")
		.default_value(std::string("all"))
		.help("cuts-type that is used in cplex; can be \"all\", \"none\", \"disj\", \"flowc\", \"flowp\", \"gomory\", \"gub\", \"liftp\", \"mcf\", \"mir\", \"zeroh\", \"clique\", \"cover\", \"impliedl\", \"impliedg\"");

	program.add_argument("--scaling")
		.default_value(std::string("off"))
		.help("CPlex scaling; can be \"off\" or \"on\"");

	program.add_argument("--model_cuts")
		.default_value(std::string("off"))
		.help("If \"on\" cuts will be added into the model as constraints; can be \"off\" or \"on\"");

	program.add_argument("--min_dim")
		.default_value(size_t(0))
		.help("minimal dimension for a problem to proceed")
		.scan<'i', size_t>();

	program.add_argument("--max_dim")
		.default_value(size_t(0))
		.help("maximal dimension for a problem to proceed; (max_dim == 0) equivalent to max_dim = +infty.")
		.scan<'i', size_t>();

	program.add_argument("--cplex_tl")
		.default_value(0.0)
		.help("time limit (sec.) for cplex MIP solving procedures")
		.scan<'g', double>();

	program.add_argument("--cuts_tl")
		.default_value(0.0)
		.help("time limit (sec.) for ILPCuts cuts generation")
		.scan<'g', double>();

	program.add_argument("--cuts_fun")
		.default_value(std::string("mir"))
		.help("manages cutting generation function; can be \"split\", \"mir\", \"pureint\" or \"hermite\"");

	program.add_argument("--trials_num")
		.default_value(size_t(5))
		.help("number of trials to choose rows for cut generation")
		.scan<'i', size_t>();

	try
	{
		program.parse_args(argc, argv);
	}
	catch (const std::exception &err)
	{
		cerr << err.what() << endl;
		cerr << program << endl;
	}

	logger log(program.get<std::string>("--output"));

	log.print_command_line_params(argc, argv);

	try
	{
		print_numeric_types_info<int_type, lint_type, rational_type, lrational_type>(log);

		check_input_params(
			log,
			program.get<double>("--tolerance"),
			program.get<double>("--err_tolerance"),
			program.get<size_t>("--cone_limit"),
			program.get<size_t>("--threads"),
			program.get<string>("--lazy"),
			program.get<int>("--node_lim"),
			program.get<string>("--cgen"),
			program.get<size_t>("--cprime"),
			program.get<string>("--cplex_cuts"),
			program.get<string>("--scaling"),
			program.get<string>("--model_cuts"),
			program.get<string>("--cuts_fun"));

		size_t threads_cuts_gen = program.get<size_t>("--threads_cuts");
		size_t threads_cplex = program.get<size_t>("--threads_cplex");
		size_t threads_global = program.get<size_t>("--threads");
		if (threads_cuts_gen == 0)
			threads_cuts_gen = threads_global;
		if (threads_cplex == 0)
			threads_cplex = threads_global;

		size_t iterations_number = program.get<size_t>("--iters");
		if (program.get<bool>("--baseline"))
			iterations_number = 0;

		run_cplex_experiment(
			log,
			program.get<string>("--input"),
			program.get<string>("--cuts_output"),
			program.get<double>("--tolerance"),
			program.get<double>("--err_tolerance"),
			program.get<size_t>("--cone_limit"),
			program.get<size_t>("--rand_strategy"),
			program.get<bool>("--sum_strat"),
			threads_cuts_gen,
			threads_cplex,
			program.get<string>("--lazy"),
			iterations_number,
			program.get<int>("--node_lim"),
			program.get<string>("--cgen"),
			program.get<size_t>("--cprime"),
			program.get<string>("--cplex_cuts"),
			program.get<string>("--scaling"),
			program.get<string>("--model_cuts"),
			program.get<size_t>("--min_dim"),
			program.get<size_t>("--max_dim"),
			program.get<double>("--cplex_tl"),
			program.get<double>("--cuts_tl"),
			program.get<size_t>("--submat_limit"),
			program.get<string>("--cuts_fun"),
			program.get<size_t>("--trials_num"));
	}
	catch (const FractionalElementException &err)
	{
		// it is not error, just can not optimize
		log.out() << err.what() << endl;
		log.out() << "ILPCuts Console: failed to optimize." << endl;
		// log.out() << program << endl;
		//  it is not error, just can not optimize
		exit(0);
	}
	catch (const std::exception &err)
	{
		log.out() << err.what() << endl;
		// log.out() << program << endl;
		std::exit(1);
	}
	catch (const IloException &err)
	{
		log.out() << err << std::endl;
		// log.out() << program << endl;
		std::exit(1);
	}

	return 0;
}