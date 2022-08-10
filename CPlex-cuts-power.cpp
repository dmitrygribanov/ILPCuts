#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

//#define ILPCUTS_PARALLEL_OFF

#include "cplex_adapter.h"
#include "conic_cuts_Eigen_Zp.h"
#include "dependencies/argparse/argparse.hpp"
#include "ILPCuts_console_utils.h"
using namespace ILPCuts;

using int_type = int64_t;
using lint_type = mpz_class;
using rational_type = double;
using lrational_type = mpq_class;

using Strategies = CutsGenStrategies<int_type, rational_type>;
using CutsGen_Base = ConicCutsGenerator_Base<int_type, rational_type, lint_type, lrational_type>;
using CutsGenAdapter = class CPlexAdapter<CutsGen_Base>;
using CutsGen_Zp = ConicCutsGenerator_Zp<int_type, rational_type, lint_type, lrational_type>;

void check_input_params(
    logger &log,
    double tolerance,
    double err_tolerance,
    size_t cone_limit_per_direction,
    size_t threads_global,
    string is_lazy_param,
    int node_lim,
    size_t cprime,
    string scaling,
    string model_cuts_switch)
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

    // if (solver_type != "adj" && solver_type != "lu" && solver_type != "zp")
    // {
    //     log.out() << "Wrong value of the CLI parameter -cg(--cgen): "
    //               << solver_type << std::endl;
    //     exit(0);
    // }

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

    // if (cplex_cut_types.check_name(cplex_cuts_param))
    // {
    //     log.out() << "Wrong value of the CLI parameter --cplex_cuts_param: "
    //               << cplex_cuts_param << std::endl;
    //     exit(0);
    // }
}

struct SolutionInfo
{
    double bestLPObjVal = 0.0;
    double bestMILPObjVal = 0.0;
    double cutoff = 0.0;
    bool is_incumbent_found = false;
};

map<string, SolutionInfo> init_solutions_map()
{
    map<string, SolutionInfo> res;

    const auto &cplex_cuts_type_names = cplex_cut_types.get_names_map();
    for (auto item : cplex_cuts_type_names)
    {
        res[item.first];
    }
    res["none"];
    res["all"];
    return res;
}

void run_cplex_experiment(
    logger &log,
    string input_path,
    string cuts_output_path,
    string comparison_result_output_path,
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
    size_t cprime,
    string scaling,
    string model_cuts_switch,
    size_t min_dim,
    size_t max_dim,
    double cplex_tl,
    double cuts_tl)
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

    // put the relaxed model into the solver
    cplex_relaxed.extract(relaxed_model);
    env.out() << "solving relaxed model ...\n";

    // run solver
    if (!cplex_relaxed.solve())
    {
        env.out() << "LP optimal solution is not found..." << endl;
        env.out() << "Solution status is " << cplex_relaxed.getStatus() << endl;
        exit(0);
    }
    double originalRelaxedObjVal = cplex_relaxed.getObjValue();

    // createing adapter object
    CutsGenAdapter cuts_gen_adapter;

    cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_Zp>(log);

    auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

    // initialization of cuts generator and setting parameters
    std::mt19937 engine;
    cuts_gen_ptr->tolerance = tolerance;
    cuts_gen_ptr->err_tolerance = err_tolerance;
    cuts_gen_ptr->cone_limit_per_direction = cone_limit_per_direction;
    cuts_gen_ptr->threads = threads_cuts_gen;
    cuts_gen_adapter.set_lazy(is_lazy_param);
    cuts_gen_ptr->rgen = &engine;
    static_cast<CutsGen_Zp *>(cuts_gen_ptr)->max_conservative_prime = cprime;

    // set strategies for cutting generation
    if (rand_directions > 0)
    {
        using RandStrategy = typename Strategies::ConvUniform<std::mt19937>;
        cuts_gen_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
    }
    if (some_direction)
        cuts_gen_ptr->direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());

    // setting parameters for solving in integer variables
    cplex_master.setParam(IloCplex::Param::Threads, threads_cplex);
    cplex_master.setParam(IloCplex::Param::MIP::Strategy::Search, CPX_MIPSEARCH_TRADITIONAL);

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

    // initialize the map, where we will store all solutions
    // corresponding to all cuts configurations
    // + "all" and "none" configurations (corresp. to all cuts and no cuts)
    auto solutions_map = init_solutions_map();

    // run cplex solver for all standard cplex cut types
    // + "all" and "none" configurations (corresp. to all cuts and no cuts)
    for (auto &solutions_iter : solutions_map)
    {
        auto &solution_info = solutions_iter.second;
        auto cut_name = solutions_iter.first;

        cplex_master.clearModel();
        cplex_master.extract(master_model);

        cplex_cut_types.set_param(cplex_master, cut_name);

        if (cplex_master.solve())
        {
            solution_info.is_incumbent_found = true;
            solution_info.bestMILPObjVal = cplex_master.getObjValue();
            solution_info.cutoff = cplex_master.getCutoff();
        }
        else
        {
            solution_info.is_incumbent_found = false;
        }
        solution_info.bestLPObjVal = cplex_master.getBestObjValue();
    }

    // construct lyambda function to add hybrid solutions to the "solutions_map"
    // cplex_param -- stores cplex cut type (+ "all" or "none")
    // solution_name -- stores hybrid configuration name
    auto add_hybrid_solution = [&](string solution_name, string cplex_param)
    {
        cplex_master.clearModel();
        cplex_cut_types.set_param(cplex_master, cplex_param);

        // creating copies of the relaxed model and constraints
        // we need this, because root_node_boosting_console_subprogram
        // write to them
        IloModel relaxed_model_copy(env);
        IloRangeArray rngs_copy(env);

        // import data to the copy of the relaxed model
        // and copy of constraints
        rngs_copy.add(rngs);
        relaxed_model_copy.add(master_model);

        // make all variables of the copy to be rational
        for (size_t i = 0; i < vars.getSize(); ++i)
        {
            relaxed_model_copy.add(IloConversion(env, vars[i], ILOFLOAT));
        }

        // solving the problem in the root node with ILPCuts cutting planes
        root_node_boosting_console_subprogram(env, cplex_master, cplex_relaxed,
                                              master_model, relaxed_model_copy,
                                              obj, rngs_copy, vars,
                                              cuts_gen_adapter,
                                              cuts_tl, iterations_number, model_cuts_switch, cuts_output_path);

        solutions_map[solution_name];

        if (cplex_master.solve())
        {
            solutions_map[solution_name].is_incumbent_found = true;
            solutions_map[solution_name].bestMILPObjVal = cplex_master.getObjValue();
            solutions_map[solution_name].cutoff = cplex_master.getCutoff();
        }
        else
        {
            solutions_map[solution_name].is_incumbent_found = false;
        }
        solutions_map[solution_name].bestLPObjVal = cplex_master.getBestObjValue();
    };

    // adding some hybrid solutions-configurations
    add_hybrid_solution("my", "none");
    add_hybrid_solution("my+all", "all");

    // printing the result into comparison_results
    env.out() << "\n";

    stringstream comparison_results;

    for (auto &solutions_iter : solutions_map)
    {
        auto &solution_info = solutions_iter.second;
        auto cut_name = solutions_iter.first;

        if (solution_info.is_incumbent_found)
        {

            comparison_results << "Cut-type = \"" << cut_name << "\": "
                               << "LP Diff. = " << originalRelaxedObjVal - solution_info.bestLPObjVal << "; "
                               << "MILP Diff. = " << originalRelaxedObjVal - solution_info.bestMILPObjVal << "; "
                               << "Cutoff Diff. = " << originalRelaxedObjVal - solution_info.cutoff << "\n"
                               << endl;
        }
        else
        {
            comparison_results << "Cut-type = \"" << cut_name << "\": No incumbent; "
                               << "LP Diff. = " << originalRelaxedObjVal - solution_info.bestLPObjVal << "; "
                               << endl;
        }
    }

    // outputing the results to all output sourses
    env.out() << comparison_results.str();

    ofstream comp_results_file(comparison_result_output_path);
    comp_results_file << comparison_results.str();
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

    program.add_argument("--comp_result")
        .default_value(std::string("comp_result.txt"))
        .help("output file with the comparison results");

    program.add_argument("--cuts_output")
        .default_value(std::string("cuts_export.txt"))
        .help("cuts export file; default is \"cuts_export.txt\"");

    // program.add_argument("--cgen")
    //     .default_value(std::string("zp"))
    //     .help("type of cutting generation algorithm; can be \"adj\", \"lu\" or \"zp\"");

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

    // program.add_argument("-b", "--baseline")
    //     .default_value(false)
    //     .help("use baseline solver only")
    //     .implicit_value(true);

    program.add_argument("--lazy")
        .default_value(string("off"))
        .help("manages CPlex lazy constraints; can be \"on\", \"off\" or \"hybrid\"");

    program.add_argument("--iters")
        .default_value(size_t(1))
        .help("number of boosting iterations in the root-node")
        .scan<'i', size_t>();

    program.add_argument("--node_lim")
        .default_value(int(1))
        .help("maximal number of nodes in CPlex B&C-algorithm (-1 = default CPlex behavor)")
        .scan<'i', int>();

    // program.add_argument("--cplex_cuts")
    //     .default_value(std::string("all"))
    //     .help("cuts-type that is used in cplex; can be \"all\", \"none\", \"disj\", \"flowc\", \"flowp\", \"gomory\", \"gub\", \"liftp\", \"mcf\", \"mir\", \"zeroh\", \"clique\", \"cover\", \"impliedl\", \"impliedg\"");

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
            program.get<size_t>("--cprime"),
            program.get<string>("--scaling"),
            program.get<string>("--model_cuts"));

        size_t threads_cuts_gen = program.get<size_t>("--threads_cuts");
        size_t threads_cplex = program.get<size_t>("--threads_cplex");
        size_t threads_global = program.get<size_t>("--threads");
        if (threads_cuts_gen == 0)
            threads_cuts_gen = threads_global;
        if (threads_cplex == 0)
            threads_cplex = threads_global;

        run_cplex_experiment(
            log,
            program.get<string>("--input"),
            program.get<string>("--cuts_output"),
            program.get<string>("--comp_result"),
            program.get<double>("--tolerance"),
            program.get<double>("--err_tolerance"),
            program.get<size_t>("--cone_limit"),
            program.get<size_t>("--rand_strategy"),
            program.get<bool>("--sum_strat"),
            threads_cuts_gen,
            threads_cplex,
            program.get<string>("--lazy"),
            program.get<size_t>("--iters"),
            program.get<int>("--node_lim"),
            program.get<size_t>("--cprime"),
            program.get<string>("--scaling"),
            program.get<string>("--model_cuts"),
            program.get<size_t>("--min_dim"),
            program.get<size_t>("--max_dim"),
            program.get<double>("--cplex_tl"),
            program.get<double>("--cuts_tl"));
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