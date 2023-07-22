#pragma once

#include <vector>

#include "numeric_module.h"

namespace ILPCuts
{
    class local_instance
    {
    public:
        enum class var_ILP_type
        {
            Int,
            Rat
        };

        is_mat A;
        // is_vec b;
        rs_vec c;
        rs_vec v;
        std::vector<var_ILP_type> ILP_var_types;

        size_t get_num_vars() const { return A.cols(); }
        size_t get_num_constr() const { return A.rows(); }
    };
}