#pragma once

#include "basic_utils.h"
#include "Eigen_solving_utilities.h"

// #include <highs/HConfig.h>
// #include <highs/util/HighsCDouble.h>
// #include <highs/util/HighsInt.h>
#include "util/HighsIntegers.h"

namespace ILPCuts
{
    namespace ModuleHiGHS
    {
        double integral_scale(const SparseVector<rational_type> &sparse_vec, double deltadown, double deltaup)
        {
            size_t n = sparse_vec.size();

            auto dense_vec = conv<SparseVector<rational_type>, vector<double>>::make(sparse_vec);

            return HighsIntegers::integralScale(dense_vec, deltadown, deltaup);
        }

        std::vector<double> integral_scale(const SparseMatrix<rational_type> &sparse_mat, double deltadown, double deltaup)
        {
            size_t rows_num = sparse_mat.rows();

            std::vector<double> result;

            for (size_t i = 0; i < rows_num; ++i)
            {
                result.push_back(integral_scale(SparseVector<rational_type>(sparse_mat.row(i)), deltadown, deltaup));
            }

            return result;
        }

        SparseMatrix<rational_type> matrix_integizer(SparseMatrix<rational_type> A, std::vector<size_t> int_col_idxs, double deltadown, double deltaup)
        {
            auto B = Utility::generate_column_submatrix_byidx<rational_type>(A, int_col_idxs);

            auto scale_vec = integral_scale(B, deltadown, deltaup);

            size_t rows_num = scale_vec.size();

            for (size_t i = 0; i < rows_num; ++i)
            {
                A.row(i) *= scale_vec[i];
            }

            return std::move(A);
        }

        SparseMatrix<rational_type> matrix_integizer(SparseMatrix<rational_type> A, std::vector<size_t> int_col_idxs, double deltadown)
        {
            return matrix_integizer(A, int_col_idxs, deltadown, deltadown);
        }
    } // namespace ModuleHiGHS
} // namespace ILPCuts