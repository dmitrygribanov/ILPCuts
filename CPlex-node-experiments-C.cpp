#include <ilcplex/ilocplex.h>
#include <ilcplex/cplexx.h>
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

int myCutCallback(CPXCENVptr env,
                  void *cbdata,
                  int wherefrom,
                  void *cbhandle,
                  int *useraction_p)
{
    int status;
    *useraction_p = CPX_CALLBACK_DEFAULT;

    // CPX_CALLBACK_MIP_CUT_LAST indicates that CPLEX is done adding cuts and the user has a last chance to add cuts
    if (wherefrom == CPX_CALLBACK_MIP_CUT_LAST)
    {
        // With this, we stop the optimization after the callback is executed
        *useraction_p = CPX_CALLBACK_FAIL;

        CPXLPptr tmpLp, copy;
        status = CPXXgetcallbacknodelp(env, cbdata, wherefrom, &tmpLp);
        copy = CPXXcloneprob(env, tmpLp, &status);
        if (status)
        {
            fprintf(stderr, "Failed to clone LP.\n");
        }

        CPXDIM precols;
        char *prectype;
        precols = CPXXgetnumcols(env, copy);
        prectype = static_cast<char *>(malloc(precols * sizeof(char)));
        status = CPXXgetcallbackctype(env, cbdata, wherefrom, prectype, 0, precols - 1);
        status = CPXXcopyctype(env, copy, prectype);
        status = CPXXwriteprob(env, copy, static_cast<char *>(cbhandle), NULL);
        if (status)
        {
            fprintf(stderr, "Failed to write presolve/cut MIP.\n");
        }
    }

    return 0;
}

string create_err_string(string prefix, CPXENVptr env, int status)
{
    char errmsg_C[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring(env, status, errmsg_C);

    return prefix + errmsg_C;
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

    status = CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    if (status)
    {
        throw runtime_error(create_err_string("Failure to turn on screen indicator...\n", env, status));
    }

    status = CPXXsetlogfilename(env, log.get_log_file_name().c_str(), "a");

    if (status)
    {
        throw runtime_error(create_err_string("Failure to set log file name into CPlex...\n", env, status));
    }

    // set global presision for long arithmetics
    mpf_set_default_prec(256);

    // create main model object
    mip = CPXXcreateprob(env, &status, "MainILP");

    if (mip == nullptr)
    {
        throw runtime_error(create_err_string("Could not create CPlex problem...\n", env, status));
    }

    status = CPXXreadcopyprob(env, mip, input_path.c_str(), nullptr);

    if (status)
    {
        throw runtime_error(create_err_string("Could not read the ILP problem...\n", env, status));
    }

    if (!baseline)
    {
        string callback_info_file_name = "reduced_model.txt";
        status = CPXXsetusercutcallbackfunc(env, myCutCallback, const_cast<char *>(callback_info_file_name.c_str()));

        if (status)
        {
            throw runtime_error(create_err_string("Could not set callback into the CPlex solver...\n", env, status));
        }

        status = CPXXsetcntparam(env, CPXPARAM_MIP_Limits_Nodes, 1);
    }

    status = CPXXsetcntparam(env, CPX_PARAM_THREADS, static_cast<int>(threads_cplex));

    status = CPXXmipopt(env, mip);

    if (status)
    {
        throw runtime_error(create_err_string("Could not optimize the ILP problem...\n", env, status));
    }

    int solstat = CPXXgetstat(env, mip);
    log.out() << "Optimization is finished!\n"
              << "Solution status is " << solstat << endl;

    // // create objective function object
    // IloObjective obj;

    // // object for variables
    // IloNumVarArray vars(env);

    // // object for inequalities
    // IloRangeArray rngs(env);

    // env.out() << "reading " << input_path << endl;

    // // import main problem from file
    // cplex_main.importModel(master_model, input_path.c_str(), obj, vars, rngs);

    // // put the original model into the solver
    // cplex_main.extract(master_model);

    // // create relaxed model object
    // IloModel relaxed_model(env);

    // // import data to the relaxed model
    // relaxed_model.add(master_model);

    // // make all variables of the relaxed model to be rational
    // for (size_t i = 0; i < vars.getSize(); ++i)
    // {
    //     relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
    // }

    // // put the relaxed model into the solver
    // cplex_relaxed.extract(relaxed_model);
    // env.out() << "solving relaxed model ...\n";

    // // setting parameters to solve relaxed model
    // // cplex_relaxed.setParam(IloCplex::Param::MIP::Interval, 1000);
    // cplex_relaxed.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

    // // run solver
    // if (!cplex_relaxed.solve())
    // {
    //     env.out() << "Optimization is terminated." << endl;
    //     env.out() << "Solution status is " << cplex_relaxed.getStatus() << endl;
    //     exit(0);
    // }
    // double originalRelaxedObjVal = cplex_relaxed.getObjValue();
    // double boostedRelaxedObjVal = 0;
    // env.out() << "solution value  is " << originalRelaxedObjVal << endl;

    // //------------------------------------------------------------------------------------------

    // // createing adapter object
    // CutsGenAdapter cuts_gen_adapter;

    // if (!baseline)
    // {

    //     cuts_gen_adapter.generator_ptr = std::make_unique<CutsGen_Zp>(log);

    //     auto cuts_gen_ptr = cuts_gen_adapter.generator_ptr.get();

    //     // initialization of cuts generator and setting parameters
    //     std::mt19937 engine;
    //     cuts_gen_ptr->rgen = &engine;
    //     cuts_gen_ptr->tolerance = tolerance;
    //     cuts_gen_ptr->err_tolerance = err_tolerance;
    //     cuts_gen_ptr->cone_limit_per_direction = cone_limit_per_direction;
    //     cuts_gen_ptr->threads = threads_cuts_gen;
    //     cuts_gen_adapter.is_lazy = CutsGenAdapter::IsLazy::on;
    //     static_cast<CutsGen_Zp *>(cuts_gen_ptr)->max_conservative_prime = cprime;

    //     // set strategies for cutting generation
    //     if (rand_directions > 0)
    //     {
    //         using RandStrategy = typename Strategies::ConvUniform<std::mt19937>;
    //         cuts_gen_ptr->direction_strategies.push_back(std::make_unique<RandStrategy>(engine, rand_directions));
    //     }
    //     if (some_direction)
    //         cuts_gen_ptr->direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());

    //     // generate cuts from the relaxed model
    //     env.out() << "constraints generation ..." << endl;
    //     cuts_gen_adapter.generate_cuts(cplex_relaxed, relaxed_model, obj, rngs, vars);

    //     // adding cuts to the original model
    //     // env.out() << "adding cuts to original model ..." << endl;
    //     // cuts_gen_adapter.add_cuts(cplex_main, vars);

    //     // setting user cuts callback
    //     // cplex_main.use(myCutsCallback(env, vars, rngs, cuts_gen_adapter));
    // }

    // //------------------------------------------------------------------------------------------

    // // setting callback
    // // int status = CPXXsetusercutcallbackfunc(env.getImpl(), cutCallbackShareRootCuts, &filename);

    // // setting parameters for solving in integer variables
    // cplex_main.setParam(IloCplex::Param::Threads, threads_cplex);
    // cplex_main.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

    // //---------------------------------------------------------------------------------------------------

    // // solve the ILP model
    // if (!cplex_main.solve())
    // {
    //     env.out() << "Optimization is terminated." << endl;
    //     env.out() << "Solution status is " << cplex_main.getStatus() << endl;
    //     exit(0);
    // }

    // env.out() << "ILP solution value  is " << cplex_main.getObjValue() << endl;
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