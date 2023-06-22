#pragma once

#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/SparseQR>

#include "basic_utils.h"

namespace ILPCuts
{
    namespace Utility
    {
        class EigenSolveException : public Exception
        {
        public:
            Eigen::ComputationInfo comp_info;

            EigenSolveException(Eigen::ComputationInfo a_comp_info)
                : comp_info(a_comp_info), Exception("ILPCuts lib error: unable to solve/decompose sparse system using Eigen.") {}

            EigenSolveException(std::string msg, Eigen::ComputationInfo a_comp_info)
                : comp_info(a_comp_info), Exception(msg) {}
        };

        template <typename ratMatrixT, typename solverT>
        ratMatrixT compute_inverse(const ratMatrixT &A, solverT &solver);

        template <typename rationalT, typename solverT>
        SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A, solverT &solver);

        template <typename rationalT>
        std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A);

        template <typename rationalT>
        rationalT compute_abs_det(const SparseMatrix<rationalT> &A);

        template <typename rationalT>
        SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b);

        template <typename rationalT, typename solverT>
        SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, solverT &solver);

        template <typename intT>
        intT vec_max2pow(const SparseVector<intT> &vec);

        template <>
        mpz_class vec_max2pow(const SparseVector<mpz_class> &vec);

        template <>
        mpz_class vec_max2pow(const SparseVector<mpz_class> &vec)
        {
            using int_type = mpz_class;
            using is_vec = SparseVector<mpz_class>;

            int_type res = 0;

            mp_bitcnt_t min_bit_index = std::numeric_limits<mp_bitcnt_t>::max();
            for (typename is_vec::InnerIterator it(vec); it; ++it)
            {
                if (it.value() != 0)
                {
                    mp_bitcnt_t bit_index = mpz_scan1(it.value().get_mpz_t(), 0);

                    // std::cout << it.index() <<", " << it.value() << std::endl;
                    min_bit_index = std::min(min_bit_index, bit_index);
                }
            }

            mpz_setbit(res.get_mpz_t(), min_bit_index);

            return std::move(res);
        }

        template <typename rationalT, typename solverT>
        SparseMatrix<rationalT> compute_inverse(solverT &solver)
        {
            using rat_mat = SparseMatrix<rationalT>;

            size_t n = solver.cols();

            rat_mat I(n, n);
            I.reserve(VectorXi::Constant(n, 2));
            for (size_t i = 0; i < n; ++i)
                I.insert(i, i) = 1;

            return std::move(solver.solve(I));
        }

        template <typename rationalT>
        SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A)
        {
            SparseLU<SparseMatrix<rationalT>, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            return compute_inverse<rationalT, decltype(decomp_LU)>(decomp_LU);
        }

        template <typename rationalT>
        std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A)
        {

            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            rat_type det = abs(decomp_LU.determinant());
            rat_mat adjA = compute_inverse<rat_type>(A, decomp_LU) * det;

            return std::make_tuple(std::move(adjA), std::move(det));
        }

        template <typename rationalT>
        rationalT compute_abs_det(const SparseMatrix<rationalT> &A)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            rat_type det = abs(decomp_LU.determinant());

            return std::move(det);
        }

        template <typename rationalT>
        size_t rank(const SparseMatrix<rationalT> &A)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseQR<rat_mat, COLAMDOrdering<int>> decomp_QR(A);

            if (decomp_QR.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_QR.info());

            return decomp_QR.rank();
        }

        template <typename rationalT>
        SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            return std::move(system_solve(b, decomp_LU));
        }

        template <typename rationalT, typename solverT>
        SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, solverT &solver)
        {
            return std::move(solver.solve(b));
        }

        template <typename T>
        SparseMatrix<T> generate_square_submatrix_byidx(const SparseMatrix<T> &A, const std::vector<size_t> &rows, const std::vector<size_t> &cols)
        {
            assert(rows.size() == cols.size());

            size_t n = rows.size();

            SparseMatrix<T> R(n, n);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < n; ++j)
                {
                    R.insert(i, j) = A.coeff(rows[i], cols[j]);
                }
            }

            R.makeCompressed();

            return R;
        }

        template <typename T>
        std::tuple<std::vector<size_t>, std::vector<size_t>> generate_some_square_submatrix(const SparseMatrix<T> &A, size_t m)
        {
            size_t n = A.cols();

            assert(n == A.rows());
            assert(m <= n);

            using rat_type = T;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            std::vector<size_t> rows(m), cols(m);

            rows = decomp_LU.rowsPermutation().indices();
            cols = decomp_LU.colsPermutation().indices();

            return std::make_tuple(rows, cols);
        }
    } // namespace Utility

} // namespace ILPCuts
