#pragma once

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "boost_tee_logging.h"

#include "cplex_adapter.h"

namespace ILPCuts
{
    using namespace std;

    template <typename intT, typename lintT, typename ratT, typename lratT>
    void print_numeric_types_info(logger &log)
    {
        log.out() << "ILPCuts numeric-types info:\n"
                  << "short-int type: " << typeid(intT).name() << "\n"
                  << "long-int type: " << typeid(lintT).name() << "\n"
                  << "short-rational type: " << typeid(ratT).name() << "\n"
                  << "long-rational type: " << typeid(lratT).name() << "\n"
                  << endl;
    }

    struct RootNodeBoostingConsoleSubprogramInfo
    {
        double total_cuts_generation_time = 0.0;
        double total_working_time = 0.0;

        double original_relaxed_objective = 0.0;
        double boosted_objective_change = 0.0;

        size_t total_iterations_num = 0;
    };

    RootNodeBoostingConsoleSubprogramInfo
    root_node_boosting_console_subprogram(IloCplex cplex_master, IloCplex cplex_relaxed,
                                          IloObjective obj, IloRangeArray orig_rngs, IloNumVarArray vars,
                                          ModuleCPlex::CPlexAdapter &cuts_gen_adapter,
                                          double cuts_tl, size_t iterations_number, string model_cuts_switch, string cuts_output_path)
    {
        IloEnv env = cplex_master.getEnv();

        cplex_relaxed.clearModel();
        cplex_master.clearModel();

        IloRangeArray rngs(env);
        rngs.add(orig_rngs);

        IloModel master_model(env);
        master_model.add(vars);
        master_model.add(rngs);
        master_model.add(obj);

        IloModel relaxed_model(env);
        relaxed_model.add(master_model);

        // make all variables of the relaxed model to be rational
        for (size_t i = 0; i < vars.getSize(); ++i)
        {
            relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
        }

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
        double boostedRelaxedObjVal = 0.0;
        double total_cut_generation_time = 0.0;
        double total_boosting_time = 0.0;

        env.out() << "solution value  is " << originalRelaxedObjVal << endl;

        if (model_cuts_switch == "off")
        {
            // put the original model into the solver
            cplex_master.extract(master_model);
        }

        if (iterations_number != 0)
        {

            env.out() << "cutting generation started..." << endl;

            // do some iterations to add cuts as constraints into the relaxed model
            // simultaniously we add cuts to the main model
            Timer boosting_timer;
            boosting_timer.start();

            for (size_t iteration_idx = 0; iteration_idx < iterations_number; ++iteration_idx)
            {
                env.out() << "root node boosting iteration " << iteration_idx << " is started..." << endl;

                // generate cuts from the relaxed model
                env.out() << "cuts generation ..." << endl;

                bool is_successed = false;
                double current_cut_generation_time = 0.0;
                std::tie(is_successed, current_cut_generation_time) = generate_cuts_with_time_limit(cuts_gen_adapter, cuts_tl, cplex_relaxed, obj, rngs, vars);

                total_cut_generation_time += current_cut_generation_time;

                if (!is_successed)
                {
                    env.out() << "\ncutting plane generation has been stopped due to the time limit (" << cuts_tl << "s),\n"
                              << "terminate..." << endl;
                    exit(0);
                }

                // checking cuts validity
                // env.out() << "initial validity checking of the generated cuts ..." << endl;
                // cuts_gen_adapter.check_pure_int_cuts_validity();

                // adding cuts to the relaxed model as constraints
                env.out() << "adding constraints to the relaxed model ..." << endl;
                IloRangeArray cuts_as_constr = cuts_gen_adapter.create_IloConstraints_from_cuts(env, vars);
                rngs.add(cuts_as_constr);
                relaxed_model.add(cuts_as_constr);

                if (model_cuts_switch == "on")
                {
                    // adding cuts to the original model as constraints
                    env.out() << "adding constraints to the original model ..." << endl;
                    master_model.add(cuts_as_constr);
                }
                else
                {
                    // adding cuts to the original model
                    env.out() << "adding cuts to original model ..." << endl;
                    cuts_gen_adapter.add_cuts(cplex_master, vars);
                }
                cuts_as_constr.end();

                // exporting cuts to a file
                if (cuts_output_path != "no")
                {
                    ofstream cuts_export_file(cuts_output_path, std::ios::out | std::ios::ate);

                    env.out() << "exporting cuts to a file ..." << endl;
                    cuts_gen_adapter.export_cuts(cuts_export_file);
                }

                // put relaxed model to the solver
                cplex_relaxed.extract(relaxed_model);
                env.out() << "solving boosted model ...\n";

                // solve relaxed model
                if (!cplex_relaxed.solve())
                {
                    env.out() << "Optimization is terminated." << endl;
                    env.out() << "Solution status is " << cplex_relaxed.getStatus() << endl;
                    exit(0);
                }
                boostedRelaxedObjVal = cplex_relaxed.getObjValue();

                env.out() << "root node boosting iteration (" << iteration_idx << ") is finished\n"
                          << "elapsed time is " << boosting_timer.finish_str()
                          << "current relaxed ILP solution value: " << boostedRelaxedObjVal << "\n"
                          << "original LP optimum: " << originalRelaxedObjVal << "\n"
                          << "boosted object value change: " << std::abs(boostedRelaxedObjVal - originalRelaxedObjVal) << "\n"
                          << "cuting planes generation time: " << current_cut_generation_time << endl;
            }

            total_boosting_time = boosting_timer.end();

            env.out() << "Cutting plane generation is finished!\n"
                      << "Number of performed iterations: " << iterations_number << "\n"
                      << "Total time (with boosting) is " << total_boosting_time << "'s\n"
                      << "original LP optimum: " << originalRelaxedObjVal << "\n"
                      << "boosted object value change: " << std::abs(boostedRelaxedObjVal - originalRelaxedObjVal) << "\n"
                      << "total cuting planes generation time: " << total_cut_generation_time << "'s" << endl;
        }

        if (model_cuts_switch == "on")
        {
            // put the original model into the solver
            cplex_master.extract(master_model);
        }

        return RootNodeBoostingConsoleSubprogramInfo{total_cut_generation_time, total_boosting_time,
                                                     originalRelaxedObjVal, std::abs(boostedRelaxedObjVal - originalRelaxedObjVal), iterations_number};
    }

    void cuts_applier_subprogram(IloCplex cplex_master,
                                 IloObjective obj, IloRangeArray rngs, IloNumVarArray vars,
                                 ModuleCPlex::CPlexAdapter &cuts_gen_adapter,
                                 string model_cuts_switch)
    {
        IloEnv env = cplex_master.getEnv();

        IloModel master_model(env);
        master_model.add(vars);
        master_model.add(rngs);
        master_model.add(obj);

        if (model_cuts_switch == "off")
        {
            cplex_master.extract(master_model);
            cuts_gen_adapter.add_cuts(cplex_master, vars);
        }
        else
        {
            IloRangeArray cuts_as_constraints = cuts_gen_adapter.create_IloConstraints_from_cuts(env, vars);
            master_model.add(cuts_as_constraints);
            cplex_master.extract(master_model);
        }
    }

    void hard_cuts_validation_subprogram(IloCplex cplex_relaxed,
                                         IloRangeArray rngs, IloNumVarArray vars,
                                         ModuleCPlex::CPlexAdapter &cuts_gen_adapter)
    {
        IloEnv env = cplex_relaxed.getEnv();

        const auto &cuts = cuts_gen_adapter.get_cuts();

        IloModel relaxed_model(env);
        relaxed_model.add(vars);
        relaxed_model.add(rngs);

        // make all variables of the relaxed model to be rational
        for (size_t i = 0; i < vars.getSize(); ++i)
        {
            relaxed_model.add(IloConversion(env, vars[i], ILOFLOAT));
        }

        for (size_t i = 0; i < cuts.size(); ++i)
        {
            env.out() << "Validating the cut with index (" << i << ")...\n";

            const auto &ray = cuts[i].ray;
            const auto &rhs = cuts[i].rhs;

            IloObjective obj = cuts_gen_adapter.apply_objective(ray, relaxed_model, vars);

            cplex_relaxed.clearModel();
            cplex_relaxed.extract(relaxed_model);

            // run solver
            if (!cplex_relaxed.solve())
            {
                env.out() << "LP optimal solution does not found...\n"
                          << "Solution status is " << cplex_relaxed.getStatus() << "\n"
                          << "Validation has failed!" << endl;
                exit(0);
            }

            double obj_val = cplex_relaxed.getObjValue();

            env.out()
                << "Max objective value is " << obj_val << "; RHS is " << rhs << "\n";
            if (static_cast<int64_t>(floor(obj_val)) != static_cast<int64_t>(rhs))
            {
                env.out() << "RHS of the cut is not correct...\n"
                          << "Abs. diff is " << abs(static_cast<double>(rhs) - obj_val) << "\n"
                          << "Validation has failed!" << endl;
                exit(0);
            }

            relaxed_model.remove(obj);
            obj.end();
        }

        env.out() << "Validation was successful!" << endl;
    }
}