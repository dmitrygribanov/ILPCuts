#pragma once

#include <vector>

namespace ILPCuts
{
    template <typename intT, typename ratT>
    class local_instance
    {
    public:
        SINGLE_NUMERIC_TYPES_DEFS(intT, ratT)

        enum class var_ILP_type
        {
            Int,
            Rat
        };

        is_mat A;
        is_vec b;
        rs_vec c;
        rs_vec v;
        std::vector<var_ILP_type> ILP_var_types;

        size_t get_num_vars() const { return A.ncols(); }
        size_t get_num_constr() const { return A.nrows(); }
    };
}