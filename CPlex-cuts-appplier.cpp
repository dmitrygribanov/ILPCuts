#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

// #define ILPCUTS_PARALLEL_OFF

#include "cplex_adapter.h"
#include "cuts_gen_cone_coverer.h"
#include "cone_coverer_Zp.h"
#include "dependencies/argparse/argparse.hpp"
#include "ILPCuts_console_utils.h"
using namespace ILPCuts;
using namespace ModuleCPlex;

using int_type = int64_t;
using lint_type = mpz_class;
using rational_type = double;
using lrational_type = mpq_class;

using ConeCoverer = ConeCoverer_Base<int_type, rational_type, lint_type, lrational_type>;
using CutsGen = ConeCoverCutsGenerator<int_type, rational_type, lint_type, lrational_type>;
using CutsGenAdapter = CPlexAdapter<CutsGen>;
using CC_Zp = ConeCoverer_Zp<int_type, rational_type, lint_type, lrational_type>;

void check_input_params(
    logger &log,
    size_t threads_global,
    string is_lazy_param,
    int node_lim,
    string scaling,
    string model_cuts_switch)
{
    if (threads_global == 0)
    {
        log.out() << "Wrong value of the CLI parameter -w(--workers, --threads): "
                  << threads_global << std::endl;
        exit(0);
    }

    if (is_lazy_param != "on" && is_lazy_param != "off" && is_lazy_param != "hybrid")
    {
        log.out() << "Wrong value of the CLI parameter --lazy: "
                  << is_lazy_param << endl;
        exit(0);
    }

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
}

void run_cplex_experiment(logger &log,
                          string input_path,
                          string cuts_import_path,
                          bool validation,
                          size_t threads,
                          string is_lazy_param,
                          int node_lim,
                          string scaling,
                          string model_cuts_switch,
                          double cplex_tl)
{
    CutsGenAdapter cuts_gen_adapter;
    cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen>(log);
    auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

    cuts_gen_ptr->cone_coverer_ptr = std::make_unique<CC_Zp>(log);
    auto cone_coverer_ptr = cuts_gen_ptr->cone_coverer_ptr.get();

    // initialization of cuts generator and setting parameters
    cuts_gen_adapter.set_lazy(is_lazy_param);

    ifstream import_cuts_file(cuts_import_path);
    cuts_gen_adapter.import_cuts(import_cuts_file);

    IloEnv env;
    env.setOut(log.out());
    env.setError(log.out());

    // set global presision for long arithmetics
    mpf_set_default_prec(256);

    // create main model object
    IloModel master_model(env);

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

    // setting parameters to solve relaxed model
    // cplex_relaxed.setParam(IloCplex::Param::MIP::Interval, 1000);
    cplex_relaxed.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    //----------------------------------------------------------
    if (validation)
    {
        hard_cuts_validation_subprogram(cplex_relaxed,
                                        rngs, vars,
                                        cuts_gen_adapter);
    }

    cuts_applier_subprogram(cplex_master,
                            obj, rngs, vars,
                            cuts_gen_adapter,
                            model_cuts_switch);
    //----------------------------------------------------------

    // setting parameters for solving in integer variables
    cplex_master.setParam(IloCplex::Param::Threads, threads);

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

    // run solver
    if (!cplex_master.solve())
    {
        env.out() << "MILP optimal solution is not found..." << endl;
        env.out() << "Solution status is " << cplex_master.getStatus() << endl;
        exit(0);
    }
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

    program.add_argument("--cuts_input")
        .default_value(std::string("cuts_export.txt"))
        .help("cuts import file; default is \"cuts_export.txt\"");

    program.add_argument("-w", "--workers", "--threads")
        .default_value(size_t(1))
        .help("global number of threads")
        .scan<'i', size_t>();

    program.add_argument("-v", "--validation")
        .default_value(false)
        .help("use hard validation procedure")
        .implicit_value(true);

    program.add_argument("--lazy")
        .default_value(string("off"))
        .help("manages CPlex lazy constraints; can be \"on\", \"off\" or \"hybrid\"");

    program.add_argument("--node_lim")
        .default_value(int(-1))
        .help("maximal number of nodes in CPlex B&C-algorithm (-1 = default CPlex behavor)")
        .scan<'i', int>();

    program.add_argument("--scaling")
        .default_value(std::string("off"))
        .help("CPlex scaling; can be \"off\" or \"on\"");

    program.add_argument("--model_cuts")
        .default_value(std::string("off"))
        .help("If \"on\" cuts will be added into the model as constraints; can be \"off\" or \"on\"");

    program.add_argument("--cplex_tl")
        .default_value(0.0)
        .help("time limit (sec.) for cplex MIP solving procedures")
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
            program.get<size_t>("--threads"),
            program.get<string>("--lazy"),
            program.get<int>("--node_lim"),
            program.get<string>("--scaling"),
            program.get<string>("--model_cuts"));

        run_cplex_experiment(
            log,
            program.get<string>("--input"),
            program.get<string>("--cuts_input"),
            program.get<bool>("--validation"),
            program.get<size_t>("--threads"),
            program.get<string>("--lazy"),
            program.get<int>("--node_lim"),
            program.get<string>("--scaling"),
            program.get<string>("--model_cuts"),
            program.get<double>("--cplex_tl"));
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