#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

//#define ILPCUTS_PARALLEL_OFF

#include <map>

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
using CutsGenAdapter = CPlexAdapter<int_type, rational_type>;
using CutsGen_Zp = ConicCutsGenerator_Zp<int_type, rational_type, lint_type, lrational_type>;
// the user-data-structure to put into CPLex C-style callback
struct MyCallbackData
{
    logger *log;
    // maximal node depth for cutting plane generation
    size_t node_max_depth;
    // maximal number of times
    // that the node is used for cutting generation
    size_t max_iters;
    CutsGenAdapter cuts_gen_adapter;

    // this data structure stores the number of times
    // that node with a given ID is used for cutting generation
    map<CPXINT, size_t> node_counters;

    bool is_node_permited_for_cuts(CPXINT nodeID, CPXINT nodeDepth)
    {
        if (node_counters.find(nodeID) != node_counters.end())
        {
            if (node_counters[nodeID] < max_iters)
            {
                ++node_counters[nodeID];
            }
            else
                return false;
        }
        else
        {
            node_counters[nodeID] = 1;
        }

        if (nodeDepth > node_max_depth)
            return false;

        return true;
    }

    logger::stream &out()
    {
        log->out() << "CPlex user callback: ";
        return log->out();
    }

    logger::stream &out(CPXINT nodeID, CPXINT nodeDepth)
    {
        log->out() << "CPlex user callback nodeID(" << nodeID << "), nodeDepth(" << nodeDepth << "):\n";
        return log->out();
    }
};

int myCutCallback(CPXCENVptr env,
                  void *cbdata,
                  int wherefrom,
                  void *cbhandle,
                  int *useraction_p)
{
    auto myData = static_cast<MyCallbackData *>(cbhandle);
    int status = 0;
    *useraction_p = CPX_CALLBACK_DEFAULT;
    CPXINT curr_node_depth = 0, curr_node_ID = 0;

    // CPX_CALLBACK_MIP_CUT_LAST indicates that CPLEX is done adding cuts and the user has a last chance to add cuts
    if (wherefrom == CPX_CALLBACK_MIP_CUT_LAST)
    {

        // reading a node depth
        status = CPXXgetcallbacknodeinfo(env, cbdata, wherefrom,
                                         0, CPX_CALLBACK_INFO_NODE_DEPTH, &curr_node_depth);
        if (status)
        {
            throw(Exception(create_err_string("CPLex CutCallback: failed to read curr. node depth...\n", env, status)));
        }

        // reading a node ID
        status = CPXXgetcallbacknodeinfo(env, cbdata, wherefrom,
                                         0, CPX_CALLBACK_INFO_NODE_SEQNUM, &curr_node_ID);
        if (status)
        {
            throw(Exception(create_err_string("CPLex CutCallback: failed to read curr. node ID...\n", env, status)));
        }

        // checking is node permited to cutting generation
        if (myData->is_node_permited_for_cuts(curr_node_ID, curr_node_depth))
        {
            *useraction_p = CPX_CALLBACK_SET;

            CPXLPptr nodeLP;
            status = CPXXgetcallbacknodelp(env, cbdata, wherefrom, &nodeLP);
            if (status)
            {
                throw(Exception(create_err_string("CPLex CutCallback: failed to get node-LP...\n", env, status)));
            }

            // copyLP = CPXXcloneprob(env, nodeLP, &status);
            // if (status)
            // {
            //     myData->out() << create_err_string("failed to clone node-LP...\n", env, status);
            // }

            // CPXDIM precols = 0;
            // precols = CPXXgetnumcols(env, nodeLP);

            // vector<char> var_types(precols, CPX_CONTINUOUS);
            //  char *prectype;
            //  prectype = static_cast<char *>(malloc(precols * sizeof(char)));
            //  status = CPXXgetcallbackctype(env, cbdata, wherefrom, prectype, 0, precols - 1);
            // status = CPXXcopyctype(env, copyLP, var_types.data());

            // status = CPXXchgprobtype(env, copyLP, CPXPROB_LP);
            // if (status)
            // {
            //     myData->out() << create_err_string("failed to change a node-LP problem type...\n", env, status);
            // }

            Timer timer;

            myData->out(curr_node_ID, curr_node_depth) << "solving of the node LP..." << endl;
            timer.start();

            status = CPXXlpopt(env, nodeLP);
            if (status)
            {
                throw(Exception(create_err_string("CPLex CutCallback: failed to optimize node-LP...\n", env, status)));
            }

            myData->out(curr_node_ID, curr_node_depth) << "solving time: " << timer.finish_str();
            // myData->out() << "callback file name ===== " << myData->callback_info_file_name.c_str() << endl;
            // status = CPXXwriteprob(env, copyLP, myData->callback_info_file_name.c_str(), NULL);
            // if (status)
            // {
            //     myData->out() << "Failed to write presolve/cut MIP." << endl;
            // }

            myData->out(curr_node_ID, curr_node_depth) << "starting ILPCuts cutting plane generation..." << endl;

            timer.start();

            myData->cuts_gen_adapter.generate_cuts(env, nodeLP);

            myData->out(curr_node_ID, curr_node_depth) << "generation time is " << timer.finish_str();

            myData->out(curr_node_ID, curr_node_depth) << "adding cuts into CPlex..." << endl;

            timer.start();

            myData->cuts_gen_adapter.add_cuts(env, cbdata, wherefrom);

            myData->out(curr_node_ID, curr_node_depth) << "adding time is " << timer.finish_str();
        }
    }

    return 0;
}

void check_input_params(
    logger &log,
    double tolerance,
    double err_tolerance,
    size_t cone_limit_per_direction,
    size_t threads_global,
    // string is_lazy_param,
    size_t iterations_number,
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

    if (iterations_number == 0)
    {
        log.out() << "Wrong value of the CLI parameter -it(--iterations): "
                  << iterations_number << std::endl;
        exit(0);
    }

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
                          size_t max_node_depth,
                          size_t iters,
                          int node_lim,
                          string scaling,
                          size_t min_dim,
                          size_t max_dim,
                          double cplex_tl)
{
    // define vars
    int status = 0;

    CPXENVptr env = nullptr;
    CPXLPptr mip = nullptr;

    //  create CPLEX envirement
    env = CPXXopenCPLEX(&status);
    if (env == nullptr)
    {
        log.out() << create_err_string("Could not open CPLEX environment...\n", env, status);
        exit(0);
    }

    // setting log-file into CPlex solver
    status = CPXXsetlogfilename(env, log.get_log_file_name().c_str(), "a");
    if (status)
    {
        log.out() << create_err_string("Failure to set log file name into CPlex...\n", env, status);
        exit(0);
    }

    // turn-on the CPlex console output
    status = CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    if (status)
    {
        log.out() << create_err_string("Failure to turn on screen indicator...\n", env, status);
        exit(0);
    }

    // setting the CPlex-env threads parameter
    status = CPXXsetcntparam(env, CPX_PARAM_THREADS, static_cast<int>(threads_cplex));
    if (status)
    {
        log.out() << create_err_string("Failed to set threads parameter for CPlex...\n", env, status);
        exit(0);
    }

    if (node_lim != -1)
    {
        status = CPXXsetcntparam(env, CPXPARAM_MIP_Limits_Nodes, node_lim);
        if (status)
        {
            log.out() << create_err_string("Failed to set node limit param. for CPlex...\n", env, status);
            exit(0);
        }
    }

    if (cplex_tl > 0.0)
    {
        status = CPXXsetcntparam(env, CPXPARAM_TimeLimit, cplex_tl);
        if (status)
        {
            log.out() << create_err_string("Failed to set time limit param. for CPlex...\n", env, status);
            exit(0);
        }
        // cplex_main.setParam(IloCplex::Param::DetTimeLimit, 0);
    }

    if (scaling == "off")
    {
        status = CPXXsetcntparam(env, CPXPARAM_MIP_SubMIP_Scale, -1);
        if (status)
        {
            log.out() << create_err_string("Failed to set SubMIP Scale param. for CPlex...\n", env, status);
            exit(0);
        }

        status = CPXXsetcntparam(env, CPXPARAM_Read_Scale, -1);
        if (status)
        {
            log.out() << create_err_string("Failed to set problem's matrix Scale param. for CPlex...\n", env, status);
            exit(0);
        }
    }

    status = CPXXsetcntparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
    status = CPXXsetcntparam(env, CPX_PARAM_DISJCUTS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_MCFCUTS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_CLIQUES, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_COVERS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_FLOWCOVERS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_FLOWPATHS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_FRACCUTS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_GUBCOVERS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_IMPLBD, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_LANDPCUTS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_LOCALIMPLBD, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_MIRCUTS, -1);
    status = CPXXsetcntparam(env, CPX_PARAM_ZEROHALFCUTS, -1);

    // set global presision for long arithmetics
    mpf_set_default_prec(256);

    // create the main MIP-problem object
    mip = CPXXcreateprob(env, &status, "MainILP");
    if (mip == nullptr)
    {
        log.out() << create_err_string("Could not create CPlex problem...\n", env, status);
        exit(0);
    }

    // reading the problem's data from file
    // nullptr - the data format is deduced authomatically
    status = CPXXreadcopyprob(env, mip, input_path.c_str(), nullptr);
    if (status)
    {
        log.out() << create_err_string("Could not read the ILP problem...\n", env, status);
        exit(0);
    }

    CPXDIM nrows = 0, ncols = 0;
    ncols = CPXXgetnumcols(env, mip);
    nrows = CPXXgetnumrows(env, mip);

    log.out() << "Number of columns in the original moddel = " << ncols << "\n"
              << "Number of rows in the original model = " << nrows << endl;

    // checking dimensional requrements for the problem
    if (ncols < min_dim)
    {
        log.out() << "\nDimension of the problem (" << ncols
                  << ") is less than min_dim (" << min_dim << ") parameter,\n"
                  << "terminate..." << endl;
        exit(0);
    }
    if (max_dim > 0 && ncols > max_dim)
    {
        log.out() << "\nDimension of the problem (" << ncols
                  << ") is greater than max_dim (" << max_dim << ") parameter,\n"
                  << "terminate..." << endl;
        exit(0);
    }
    // //------------------------------------------------------------------------------------------

    // creating data for my cut-callback
    MyCallbackData myData = {&log, max_node_depth, iters};

    if (!baseline)
    {
        // registering my cut-callback
        status = CPXXsetusercutcallbackfunc(env, myCutCallback, &myData);
        if (status)
        {
            log.out() << create_err_string("Could not set callback into the CPlex solver...\n", env, status);
            exit(0);
        }

        auto &cuts_gen_adapter = myData.cuts_gen_adapter;

        // creating ILPCuts generator of zp-type
        cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_Zp>(log);
        auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

        // initialization of ILPCuts generator and setting its parameters
        std::mt19937 engine;
        cuts_gen_ptr->rgen = &engine;
        cuts_gen_ptr->tolerance = tolerance;
        cuts_gen_ptr->err_tolerance = err_tolerance;
        cuts_gen_ptr->cone_limit_per_direction = cone_limit_per_direction;
        cuts_gen_ptr->threads = threads_cuts_gen;
        cuts_gen_adapter.is_lazy = CutsGenAdapter::IsLazy::off;
        static_cast<CutsGen_Zp *>(cuts_gen_ptr)->max_conservative_prime = cprime;

        // set strategies for cutting generation
        if (rand_directions > 0)
        {
            using RandStrategy = typename Strategies::ConvUniform<std::mt19937>;
            cuts_gen_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
        }
        if (some_direction)
            cuts_gen_ptr->direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());

        // status = CPXXsetcntparam(env, CPXPARAM_MIP_Limits_Nodes, 1);
    }

    // starting the optimization process
    status = CPXXmipopt(env, mip);
    if (status)
    {
        log.out() << create_err_string("Could not optimize the ILP problem...\n", env, status)
                  << endl;
    }
    else
    {
        int solstat = CPXXgetstat(env, mip);
        log.out() << "Optimization is finished!\n"
                  << "Solution status is " << solstat << endl;

        double MIP_opt_val = 0.0;
        status = CPXXgetobjval(env, mip, &MIP_opt_val);
        if (status)
        {
            log.out() << create_err_string("Could not extract the optimal value of the ILP objective...\n", env, status)
                      << endl;
        }

        log.out() << "Optimal value of the objective: " << MIP_opt_val
                  << endl;
    }
}

int main(int argc, char **argv)
{
    argparse::ArgumentParser program("CPlexNodeExpsC");

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

    program.add_argument("-d", "--max_node_depth")
        .default_value(size_t(0))
        .help("max node depth for ILPCuts-cuts generation")
        .scan<'i', size_t>();

    // program.add_argument("--lazy")
    //     .default_value(string("off"))
    //     .help("manages CPlex lazy constraints; can be \"on\", \"off\" or \"hybrid\"");

    program.add_argument("--iters")
        .default_value(size_t(1))
        .help("maximal number of iterations of cutting generation in the root-node")
        .scan<'i', size_t>();

    program.add_argument("--node_lim")
        .default_value(int(-1))
        .help("maximal number of nodes in CPlex B&C-algorithm (-1 = default CPlex behavor)")
        .scan<'i', int>();

    program.add_argument("--scaling")
        .default_value(std::string("off"))
        .help("CPlex scaling; can be \"off\" or \"on\"");

    // program.add_argument("--cplex_cuts")
    //     .default_value(std::string("all"))
    //     .help("cuts-type that is used in cplex; can be \"all\", \"none\", \"disj\", \"flowc\", \"flowp\", \"gomory\", \"gub\", \"liftp\", \"mcf\", \"mir\", \"zeroh\", \"clique\", \"cover\", \"impliedl\", \"impliedg\"");

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
        .help("time limit for cplex MIP solving procedures")
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
            program.get<size_t>("--iters"),
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
            program.get<size_t>("--max_node_depth"),
            program.get<size_t>("--iters"),
            program.get<int>("--node_lim"),
            program.get<string>("--scaling"),
            program.get<size_t>("--min_dim"),
            program.get<size_t>("--max_dim"),
            program.get<double>("--cplex_tl"));
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