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
using CutsGenAdapter = CPlexAdapter<CutsGen_Base>;
using CutsGen_Zp = ConicCutsGenerator_Zp<int_type, rational_type, lint_type, lrational_type>;

// ILOUSERCUTCALLBACK4(myCutsCallback, IloCplex, cplex, IloObjective, obj, IloNumVarArray, vars, IloRangeArray, rangs)

ILOUSERCUTCALLBACK3(myCutsCallback, IloNumVarArray, vars, IloRangeArray, rngs, CutsGenAdapter &, adapter)
{
    IloEnv env = getEnv();
    IloModel mod = getModel();

    if (getCurrentNodeDepth() == 0)
    {
        IloRangeArray ilo_cuts(env);

        // cout << "----------------" << endl;
        if (!adapter.get_cuts().empty())
        {
            adapter.transform_cuts_to_ilo_constraints(env, vars, ilo_cuts);

            for (size_t i = 0; i < ilo_cuts.getSize(); ++i)
            {
                // cout << "------" << i << "---------" << endl;
                addLocal(ilo_cuts[i]);
            }

            adapter.clear_cuts();

            ilo_cuts.end();
        }
    }

    IloIterator<IloExpr> citer = IloIterator<IloExpr>(env);
    IloModel::Iterator iter = IloModel::Iterator(mod);
    IloConstraintArray carr(env);
    IloNumVarArray varr(env);
    IloExprArray earr(env);
    size_t cc = 0, vc = 0, ec = 0;
    while (iter.ok())
    {
        if ((*iter).isConstraint())
        {
            carr.add((*iter).asConstraint());
            ++cc;
        }
        else if ((*iter).isVariable())
        {
            varr.add((*iter).asVariable());
            ++vc;
        }
        else if ((*iter).isNumExpr())
        {
            earr.add((*iter).asNumExpr());
            ++ec;
        }
        else if ((*iter).isIntExpr())
        {
            earr.add((*iter).asIntExpr());
            ++ec;
        }
        ++iter;
    }

    env.out() << "=====variables: " << vc << endl;
    // for (size_t i = 0; i < varr.getSize(); ++i)
    // {
    //     env.out() << varr[i].getName() << ", ";
    // }
    // env.out() << endl;

    env.out() << "=====constraints: " << cc << endl;
    // for (size_t i = 0; i < carr.getSize(); ++i)
    // {
    //     env.out() << carr[i].getName() << ", ";
    // }
    // env.out() << endl;

    // env.out() << "=====constraints 2.0==========" << endl;
    // size_t expr_count = 0;
    // while (citer.ok())
    // {
    //     if ((*citer).isNumExpr())
    //     {
    //         env.out() << (*citer).asNumExpr() << ", ";
    //         ++expr_count;
    //     }
    //     else if ((*citer).isIntExpr())
    //     {
    //         env.out() << (*citer).asIntExpr() << ", ";
    //         ++expr_count;
    //     }

    //     ++citer;
    // }

    // env.out()<< expr_count << '\n';
    /*
    env.out() << "=====exprs==========" << endl;
    for (size_t i = 0; i < earr.getSize(); ++i) {
        env.out() << earr[i].getName() << ", ";
    }

    */

    IloNumArray slacks_arr(env);
    IloNumArray vars_arr(env);
    IloNumArray UBs_arr(env);
    IloNumArray LBs_arr(env);
    getSlacks(slacks_arr, rngs);
    getValues(vars_arr, vars);
    getUBs(UBs_arr, vars);
    getLBs(LBs_arr, vars);

    size_t base_num = 0;
    for (size_t i = 0; i < slacks_arr.getSize(); i++)
    {
        if (slacks_arr[i] == 0)
            ++base_num;
    }
    for (size_t i = 0; i < vars.getSize(); i++)
    {
        if (vars_arr[i] == UBs_arr[i] || vars_arr[i] == UBs_arr[i])
            ++base_num;
    }
    env.out() << "vars num = " << vars.getSize()
              << "; slacks num = " << slacks_arr.getSize()
              << "; base num = " << base_num << endl;

    vars_arr.end();
    UBs_arr.end();
    LBs_arr.end();
    slacks_arr.end();
    carr.end();
    varr.end();
    earr.end();
}

void check_input_params(
    logger &log,
    double tolerance,
    double err_tolerance,
    size_t cone_limit_per_direction,
    size_t threads_global,
    // string is_lazy_param,
    // size_t iterations_number,
    int node_lim,
    // std::string solver_type,
    size_t cprime,
    // string cplex_cuts_param,
    string scaling
    // string model_cuts_switch
)
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

    // if (is_lazy_param != "on" && is_lazy_param != "off" && is_lazy_param != "hybrid")
    // {
    //     log.out() << "Wrong value of the CLI parameter --lazy: "
    //               << is_lazy_param << endl;
    //     exit(0);
    // }

    // if (iterations_number == 0)
    // {
    //     log.out() << "Wrong value of the CLI parameter -it(--iterations): "
    //               << iterations_number << std::endl;
    //     exit(0);
    // }

    // if (model_cuts_switch != "on" && model_cuts_switch != "off")
    // {
    //     log.out() << "Wrong value of the CLI parameter --model_cuts: "
    //               << model_cuts_switch << std::endl;
    //     exit(0);
    // }

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

    // if (cplex_cuts_param != "all" && cplex_cuts_param != "none" && cplex_cut_types.find(cplex_cuts_param) == cplex_cut_types.end())
    // {
    //     log.out() << "Wrong value of the CLI parameter --cplex_cuts_param: "
    //               << cplex_cuts_param << std::endl;
    //     exit(0);
    // }
}

void run_cplex_experiment(logger &log,
                          string input_path,
                          double tolerance,
                          double err_tolerance,
                          size_t cone_limit_per_direction,
                          size_t rand_directions,
                          bool some_direction,
                          size_t threads_cuts_gen,
                          size_t threads_cplex,
                          bool baseline,
                          size_t cprime,
                          int node_lim,
                          string scaling)
{
    // create CPLEX envirement
    IloEnv env;
    env.setOut(log.out());
    env.setError(log.out());

    // set global presision for long arithmetics
    mpf_set_default_prec(256);

    // create main model object
    IloModel master_model(env);

    // create solver object
    IloCplex cplex_main(env);

    // solver object to solve relaxed model
    IloCplex cplex_relaxed(env);

    // create objective function object
    IloObjective obj;

    // object for variables
    IloNumVarArray vars(env);

    // object for inequalities
    IloRangeArray rngs(env);

    env.out() << "reading " << input_path << endl;

    // import main problem from file
    cplex_main.importModel(master_model, input_path.c_str(), obj, vars, rngs);

    // put the original model into the solver
    cplex_main.extract(master_model);

    // create relaxed model object
    IloModel relaxed_model(env);

    // import data to the relaxed model
    relaxed_model.add(master_model);

    // make all variables of the relaxed model to be rational
    for (size_t i = 0; i < vars.getSize(); ++i)
    {
        relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
    }

    // put the relaxed model into the solver
    cplex_relaxed.extract(relaxed_model);
    env.out() << "solving relaxed model ...\n";

    // setting parameters to solve relaxed model
    // cplex_relaxed.setParam(IloCplex::Param::MIP::Interval, 1000);
    cplex_relaxed.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    // run solver
    if (!cplex_relaxed.solve())
    {
        env.out() << "Optimization is terminated." << endl;
        env.out() << "Solution status is " << cplex_relaxed.getStatus() << endl;
        exit(0);
    }
    double originalRelaxedObjVal = cplex_relaxed.getObjValue();
    double boostedRelaxedObjVal = 0;
    env.out() << "solution value  is " << originalRelaxedObjVal << endl;

    //------------------------------------------------------------------------------------------

    // createing adapter object
    CutsGenAdapter cuts_gen_adapter;

    if (!baseline)
    {
        cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_Zp>(log);

        auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

        // initialization of cuts generator and setting parameters
        std::mt19937 engine;
        cuts_gen_ptr->rgen = &engine;
        cuts_gen_ptr->tolerance = tolerance;
        cuts_gen_ptr->err_tolerance = err_tolerance;
        cuts_gen_ptr->cone_limit_per_direction = cone_limit_per_direction;
        cuts_gen_ptr->threads = threads_cuts_gen;
        cuts_gen_adapter.is_lazy = CutsGenAdapter::IsLazy::on;
        static_cast<CutsGen_Zp *>(cuts_gen_ptr)->max_conservative_prime = cprime;

        // set strategies for cutting generation
        if (rand_directions > 0)
        {
            using RandStrategy = typename Strategies::ConvUniform<std::mt19937>;
            cuts_gen_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
        }
        if (some_direction)
            cuts_gen_ptr->direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());

        // generate cuts from the relaxed model
        env.out() << "constraints generation ..." << endl;
        cuts_gen_adapter.generate_cuts(cplex_relaxed, relaxed_model, obj, rngs, vars);

        // adding cuts to the original model
        // env.out() << "adding cuts to original model ..." << endl;
        // cuts_gen_adapter.add_cuts(cplex_main, vars);

        // setting user cuts callback
        cplex_main.use(myCutsCallback(env, vars, rngs, cuts_gen_adapter));
    }

    //------------------------------------------------------------------------------------------

    // setting parameters for solving in integer variables
    cplex_main.setParam(IloCplex::Param::Threads, threads_cplex);
    cplex_main.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

    if (scaling != "on" && scaling != "off")
    {
        log.out() << "Wrong value of the CLI parameter --scaling: " << scaling << std::endl;
        exit(0);
    }

    if (scaling == "off")
    {
        cplex_main.setParam(IloCplex::Param::Read::Scale, -1);
        // cplex_main.setParam(IloCplex::Param::MIP::SubMIP::Scale, -1);
    }

    if (node_lim < -1)
    {
        env.error() << "Wrong value of the CLI parameter --node_lim: " << node_lim << std::endl;
        exit(0);
    }
    if (node_lim != -1)
        cplex_main.setParam(IloCplex::Param::MIP::Limits::Nodes, node_lim);

    //---------------------------------------------------------------------------------------------------

    // solve the ILP model
    if (!cplex_main.solve())
    {
        env.out() << "Optimization is terminated." << endl;
        env.out() << "Solution status is " << cplex_main.getStatus() << endl;
        exit(0);
    }

    env.out() << "ILP solution value  is " << cplex_main.getObjValue() << endl;
}

int main(int argc, char **argv)
{
    argparse::ArgumentParser program("CPlexNodeExpsCPP");

    program.add_argument("-i", "--input")
        .default_value(std::string("-"))
        .required()
        .help("specify the input file.");

    program.add_argument("-o", "--output")
        .default_value(std::string("log.txt"))
        .help("output log file; default is \"log.txt\"");

    // program.add_argument("--cuts_output")
    //     .default_value(std::string("cuts_export.txt"))
    //     .help("cuts export file; default is \"cuts_export.txt\"");

    // program.add_argument("--cgen")
    //     .default_value(std::string("adj"))
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

    program.add_argument("-b", "--baseline")
        .default_value(false)
        .help("use baseline solver only")
        .implicit_value(true);

    // program.add_argument("--lazy")
    //     .default_value(string("off"))
    //     .help("manages CPlex lazy constraints; can be \"on\", \"off\" or \"hybrid\"");

    // program.add_argument("--iters")
    //     .default_value(size_t(1))
    //     .help("number of boosting iterations in the root-node")
    //     .scan<'i', size_t>();

    program.add_argument("--node_lim")
        .default_value(int(-1))
        .help("maximal number of nodes in CPlex B&C-algorithm (-1 = default CPlex behavor)")
        .scan<'i', int>();

    // program.add_argument("--cplex_cuts")
    //     .default_value(std::string("all"))
    //     .help("cuts-type that is used in cplex; can be \"all\", \"none\", \"disj\", \"flowc\", \"flowp\", \"gomory\", \"gub\", \"liftp\", \"mcf\", \"mir\", \"zeroh\", \"clique\", \"cover\", \"impliedl\", \"impliedg\"");

    program.add_argument("--scaling")
        .default_value(std::string("off"))
        .help("CPlex scaling; can be \"off\" or \"on\"");

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

        check_input_params(log,
                           program.get<double>("--tolerance"),
                           program.get<double>("--err_tolerance"),
                           program.get<size_t>("--cone_limit"),
                           program.get<size_t>("--threads"),
                           program.get<int>("--node_lim"),
                           program.get<size_t>("--cprime"),
                           program.get<string>("--scaling"));

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
            program.get<double>("--tolerance"),
            program.get<double>("--err_tolerance"),
            program.get<size_t>("--cone_limit"),
            program.get<size_t>("--rand_strategy"),
            program.get<bool>("--sum_strat"),
            threads_cuts_gen,
            threads_cplex,
            program.get<bool>("--baseline"),
            program.get<size_t>("--cprime"),
            program.get<int>("--node_lim"),
            program.get<string>("--scaling"));
    }
    catch (const FractionalElementException &err)
    {
        // it is not error, just can not optimize
        log.out() << err.what() << endl;
        log.out() << "ILPCuts Console: failed to optimize." << endl;

        //  it is not error, just can not optimize
        exit(0);
    }
    catch (const std::exception &err)
    {
        log.out() << err.what() << endl;
        std::exit(1);
    }
    catch (const IloException &err)
    {
        log.out() << err << std::endl;
        std::exit(1);
    }

    return 0;
}