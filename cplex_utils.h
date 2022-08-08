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

    class CPlexCutTypes
    {
    private:
        map<string, IloCplex::IntParam> _names_map;

    public:
        CPlexCutTypes()
        {
            _names_map["clique"] = IloCplex::Param::MIP::Cuts::Cliques;
            _names_map["cover"] = IloCplex::Param::MIP::Cuts::Covers;
            _names_map["gomory"] = IloCplex::Param::MIP::Cuts::Gomory;
            _names_map["disj"] = IloCplex::Param::MIP::Cuts::Disjunctive;
            _names_map["flowc"] = IloCplex::Param::MIP::Cuts::FlowCovers;
            _names_map["flowp"] = IloCplex::Param::MIP::Cuts::PathCut;
            _names_map["gub"] = IloCplex::Param::MIP::Cuts::GUBCovers;
            _names_map["liftp"] = IloCplex::Param::MIP::Cuts::LiftProj;
            _names_map["mcf"] = IloCplex::Param::MIP::Cuts::MCFCut;
            _names_map["mir"] = IloCplex::Param::MIP::Cuts::MIRCut;
            _names_map["zeroh"] = IloCplex::Param::MIP::Cuts::ZeroHalfCut;
            _names_map["impliedl"] = IloCplex::Param::MIP::Cuts::LocalImplied;
            _names_map["impliedg"] = IloCplex::Param::MIP::Cuts::Implied;
        }

        const map<string, IloCplex::IntParam> &get_names_map() const
        {
            return _names_map;
        }

        bool check_name(string name) const
        {
            return name != "all" && name != "none" && _names_map.find(name) == _names_map.end();
        }

        void set_param(IloCplex cplex, string name) const
        {
            // if all-type is used, we run defoult cplex solver
            if (name != "all")
            {
                // in the opposite case, if type == none, we turn off all cplex cuts
                if (name == "none")
                {
                    for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                        cplex.setParam(types_iter->second, -1);
                }
                // if type != none, we turn off all cuts except (cplex_cuts_param)-variable
                else
                {
                    for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                        if (types_iter->first != name)
                            cplex.setParam(types_iter->second, -1);
                }
            }
        }
    };

    const CPlexCutTypes cplex_cut_types;

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
