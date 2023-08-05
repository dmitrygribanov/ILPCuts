#pragma once

#include "basic_utils.h"
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplexx.h>
#include <map>

namespace ILPCuts
{
    namespace ModuleCPlex
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
                return name != "all(force)" && name != "all" && name != "none" && _names_map.find(name) == _names_map.end();
            }

            void restore_params(IloCplex cplex) const
            {
                for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                    cplex.setParam(types_iter->second, 0);
            }

            void set_param(IloCplex cplex, string name) const
            {
                // restore params to default if they were changed before
                restore_params(cplex);

                // if all-type is used, we run defoult cplex solver
                if (name != "all")
                {
                    // in the opposite case, if type == none, we turn off all cplex cuts
                    if (name == "none")
                    {
                        for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                            cplex.setParam(types_iter->second, -1);
                    }
                    else if (name == "all(force)")
                    {
                        for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                            cplex.setParam(types_iter->second, 2);
                    }
                    // if type != none and all(force), we turn off all cuts except (cplex_cuts_param)-variable
                    // the last is set to agressive(=2) switch
                    else
                    {
                        for (auto types_iter = _names_map.begin(); types_iter != _names_map.end(); ++types_iter)
                            if (types_iter->first != name)
                                cplex.setParam(types_iter->second, -1);
                            else
                                cplex.setParam(types_iter->second, 2);
                    }
                }
            }
        };

        const CPlexCutTypes cplex_cut_types;

        class FractionalElementException : public Exception
        {
        public:
            FractionalElementException()
                : Exception("ILPCuts lib error: fractional element in A-matrix has been found.") {}

            FractionalElementException(size_t base_row_idx, size_t col_idx, double val)
                : Exception(std::string("ILPCuts lib error: fracttional elemen in A has located;\n") +
                            "ConstrNum = " + std::to_string(base_row_idx) + "\n" +
                            "VarNum = " + std::to_string(col_idx) + "\n" +
                            "value = " + std::to_string(val))
            {
            }
        };

        class SemiVariableException : public Exception
        {
        public:
            SemiVariableException()
                : Exception("ILPCuts lib error: the problem contains semi-integer/continuous variables.")
            {
            }

            SemiVariableException(size_t varID, char type)
                : Exception("ILPCuts lib error: the problem contains semi-integer/continuous variables, varID = " + std::to_string(varID) + ", type = " + type)
            {
            }
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

        string print_IloRange(IloRange rng)
        {
            stringstream res;

            res << rng.getLB() << " <= " << rng.getUB() << "\n";
            IloExpr::LinearIterator riter = rng.getLinearIterator();
            while (riter.ok())
            {
                res << riter.getCoef() << "*" << riter.getVar() << " ";
                ++riter;
            }

            return res.str();
        }

        string print_IloRangeArray(IloRangeArray rngs)
        {
            stringstream res;

            for (size_t i = 0; i < rngs.getSize(); ++i)
            {
                IloRange rng = rngs[i];
                res << "range(" << i << ") : " << print_IloRange(rng) << "\n";
            }

            return res.str();
        }
    } // namespace ModuleCPlex

} // namespace ILPCuts
