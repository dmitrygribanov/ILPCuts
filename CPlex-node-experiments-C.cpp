#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

//#define ILPCUTS_PARALLEL_OFF

#include "cplex_adapter.h"
#include "conic_cuts_Eigen_Zp.h"
#include "dependencies/argparse/argparse.hpp"
using namespace ILPCuts;

using int_type = int64_t;
using lint_type = mpz_class;
using rational_type = double;
using lrational_type = mpq_class;

using Strategies = CutsGenStrategies<int_type, rational_type>;
using CutsGen_Base = ConicCutsGenerator_Base<int_type, rational_type, lint_type, lrational_type>;
using CutsGenAdapter = CPlexAdapter<CutsGen_Base>;
using CutsGen_Zp = ConicCutsGenerator_Zp<int_type, rational_type, lint_type, lrational_type>;

string create_err_string(string prefix, CPXCENVptr env, int status)
{
    char errmsg_C[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring(env, status, errmsg_C);

    return prefix + errmsg_C;
}
struct MyCallbackData
{
    string callback_info_file_name;
    logger *log;
    CutsGenAdapter cuts_gen_adapter;

    logger::stream &out()
    {
        log->out() << "CPlex user callback: ";
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
    int status;
    *useraction_p = CPX_CALLBACK_DEFAULT;

    // CPX_CALLBACK_MIP_CUT_LAST indicates that CPLEX is done adding cuts and the user has a last chance to add cuts
    if (wherefrom == CPX_CALLBACK_MIP_CUT_LAST)
    {
        // With this, we stop the optimization after the callback is executed
        *useraction_p = CPX_CALLBACK_FAIL;

        CPXLPptr nodeLP, copyLP;
        status = CPXXgetcallbacknodelp(env, cbdata, wherefrom, &nodeLP);
        if (status)
        {
            myData->out() << create_err_string("failed to get node-LP...\n", env, status);
        }

        copyLP = CPXXcloneprob(env, nodeLP, &status);
        if (status)
        {
            myData->out() << create_err_string("failed to clone node-LP...\n", env, status);
        }

        CPXDIM precols;
        precols = CPXXgetnumcols(env, copyLP);

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

        status = CPXXlpopt(env, copyLP);
        if (status)
        {
            myData->out() << create_err_string("failed to optimize node-LP...\n", env, status);
        }
        // myData->out() << "callback file name ===== " << myData->callback_info_file_name.c_str() << endl;
        // status = CPXXwriteprob(env, copyLP, myData->callback_info_file_name.c_str(), NULL);
        // if (status)
        // {
        //     myData->out() << "Failed to write presolve/cut MIP." << endl;
        // }
    }

    return 0;
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
                          size_t cprime)
{
    // define vars
    int status = 0;

    CPXENVptr env = nullptr;
    CPXLPptr mip = nullptr;

    //  create CPLEX envirement
    env = CPXXopenCPLEX(&status);
    if (env == nullptr)
    {
        throw runtime_error(create_err_string("Could not open CPLEX environment...\n", env, status));
    }

    // turn-on the CPlex console output
    status = CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    if (status)
    {
        throw runtime_error(create_err_string("Failure to turn on screen indicator...\n", env, status));
    }

    // setting log-file into CPlex solver
    status = CPXXsetlogfilename(env, log.get_log_file_name().c_str(), "a");
    if (status)
    {
        throw runtime_error(create_err_string("Failure to set log file name into CPlex...\n", env, status));
    }

    // set global presision for long arithmetics
    mpf_set_default_prec(256);

    // create the main MIP-problem object
    mip = CPXXcreateprob(env, &status, "MainILP");
    if (mip == nullptr)
    {
        throw runtime_error(create_err_string("Could not create CPlex problem...\n", env, status));
    }

    // reading the problem's data from file
    // nullptr - the data format is deduced authomatically
    status = CPXXreadcopyprob(env, mip, input_path.c_str(), nullptr);
    if (status)
    {
        throw runtime_error(create_err_string("Could not read the ILP problem...\n", env, status));
    }

    log.out() << "Number of columns in the original moddel = " << CPXXgetnumcols(env, mip) << "\n"
              << "Number of rows in the original model = " << CPXXgetnumrows(env, mip) << endl;
    // //------------------------------------------------------------------------------------------

    // creating data for my cut-callback
    MyCallbackData myData = {"reduced_model.mps", &log};

    if (!baseline)
    {
        // registering my cut-callback
        status = CPXXsetusercutcallbackfunc(env, myCutCallback, &myData);
        if (status)
        {
            throw runtime_error(create_err_string("Could not set callback into the CPlex solver...\n", env, status));
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

    // setting the CPlex MIP-solver parameters
    status = CPXXsetcntparam(env, CPX_PARAM_THREADS, static_cast<int>(threads_cplex));
    {
        throw runtime_error(create_err_string("Failed to set threads parameter for CPlex...\n", env, status));
    }

    // starting the optimization process
    status = CPXXmipopt(env, mip);
    if (status)
    {
        throw runtime_error(create_err_string("Could not optimize the ILP problem...\n", env, status));
    }

    int solstat = CPXXgetstat(env, mip);
    log.out() << "Optimization is finished!\n"
              << "Solution status is " << solstat << endl;
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
        .default_value(size_t(101))
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

    // program.add_argument("--cplex_cuts")
    //     .default_value(std::string("all"))
    //     .help("cuts-type that is used in cplex; can be \"all\", \"none\", \"disj\", \"flowc\", \"flowp\", \"gomory\", \"gub\", \"liftp\", \"mcf\", \"mir\", \"zeroh\", \"clique\", \"cover\", \"impliedl\", \"impliedg\"");

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

    try
    {
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
            program.get<size_t>("--cprime"));
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