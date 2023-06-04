#pragma once

#include <vector>

namespace ILPCuts
{
    template <typename intT, typename ratT>
    class LocalInstance
    {
    public:
        using rational_type = ratT;
        using int_type = intT;

        using is_mat = SparseMatrix<int_type>;
        using rs_mat = SparseMatrix<rational_type>;

        using is_vec = SparseVector<int_type>;
        using rs_vec = SparseVector<rational_type>;

        enum class var_ILP_type
        {
            Int,
            Rat
        };

        is_mat A;
        is_vec b;
        rs_vec c;
        std::vector<var_ILP_type> ILP_var_types;

        size_t get_num_vars() const { return A.ncols(); }
        size_t get_num_constr() const { return A.nrows(); }
    };
}