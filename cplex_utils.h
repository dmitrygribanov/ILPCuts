#pragma once

#include "basic_utils.h"
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplexx.h>
#include <map>

namespace ILPCuts
{
    string create_err_string(string prefix, CPXCENVptr env, int status, string suff = "")
    {
        stringstream res_msg;
        res_msg << prefix;

        char buff_C[CPXMESSAGEBUFSIZE];
        const char *msg_C = CPXXgeterrorstring(env, status, buff_C);

        if (msg_C)
            res_msg << msg_C;
        else
            res_msg << "CPlex Error " << status << ": unknown error code.\n";

        res_msg << suff;

        return res_msg.str();
    }

    map<string, IloCplex::IntParam> init_cplex_cut_type_names_map()
    {
        map<string, IloCplex::IntParam> res;

        res["clique"] = IloCplex::Param::MIP::Cuts::Cliques;
        res["cover"] = IloCplex::Param::MIP::Cuts::Covers;
        res["gomory"] = IloCplex::Param::MIP::Cuts::Gomory;
        res["disj"] = IloCplex::Param::MIP::Cuts::Disjunctive;
        res["flowc"] = IloCplex::Param::MIP::Cuts::FlowCovers;
        res["flowp"] = IloCplex::Param::MIP::Cuts::PathCut;
        res["gub"] = IloCplex::Param::MIP::Cuts::GUBCovers;
        res["liftp"] = IloCplex::Param::MIP::Cuts::LiftProj;
        res["mcf"] = IloCplex::Param::MIP::Cuts::MCFCut;
        res["mir"] = IloCplex::Param::MIP::Cuts::MIRCut;
        res["zeroh"] = IloCplex::Param::MIP::Cuts::ZeroHalfCut;
        res["impliedl"] = IloCplex::Param::MIP::Cuts::LocalImplied;
        res["impliedg"] = IloCplex::Param::MIP::Cuts::Implied;

        return res;
    }
    class FractionalElementException : public Exception
    {
    public:
        FractionalElementException()
            : Exception("ILPCuts lib error: fractional element in A-matrix has been found") {}

        FractionalElementException(std::string msg)
            : Exception(msg) {}
    };

    IloNumVarArray extract_vars(IloEnv env, IloModel model)
    {
        IloModel::Iterator model_iter = IloModel::Iterator(model);

        IloNumVarArray res(env);

        while (model_iter.ok())
        {
            if ((*model_iter).isVariable())
            {
                res.add((*model_iter).asVariable());
            }

            ++model_iter;
        }

        return res;
    }

    IloObjective extract_objective(IloEnv env, IloModel model)
    {
        IloModel::Iterator model_iter = IloModel::Iterator(model);

        IloObjective res;

        size_t objective_counter = 0;

        while (model_iter.ok())
        {
            if ((*model_iter).isObjective())
            {
                if (objective_counter > 0)
                {
                    throw Exception("ILPCuts error: multople CPlex objective detected...");
                }
                else
                {
                    ++objective_counter;
                    res = (*model_iter).asObjective();
                }
            }

            ++model_iter;
        }

        return res;
    }

    IloConstraintArray extract_constraints(IloEnv env, IloModel model)
    {
        IloModel::Iterator model_iter = IloModel::Iterator(model);

        IloConstraintArray res(env);

        while (model_iter.ok())
        {
            if ((*model_iter).isConstraint())
            {
                res.add((*model_iter).asConstraint());
            }

            ++model_iter;
        }

        return res;
    }
} // namespace ILPCuts
