#pragma once

#include "base_module.h"

namespace ILPCuts
{
    namespace Utility
    {
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
            /*while (
                (vec_copy.unaryExpr(std::ref(mod2_unaryop)) == 0).all()
                ) {
                res *= 2;
                vec_copy /= 2;
            }*/

            mpz_setbit(res.get_mpz_t(), min_bit_index);

            return std::move(res);
        }

        template <typename rationalT, typename solverT>
        SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A, solverT &solver)
        {
            using rat_mat = SparseMatrix<rationalT>;

            size_t n = A.cols();

            rat_mat I(n, n);
            I.reserve(VectorXi::Constant(n, 2));
            for (size_t i = 0; i < n; ++i)
                I.insert(i, i) = 1;

            rat_mat res(n, n);
            /*for (size_t i = 0; i < n; ++i) {
                res.col(i) = solver.solve(I.col(i));
            }*/
            res = std::move(solver.solve(I));

            return std::move(res);
        }

        template <typename rationalT>
        std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A)
        {

            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

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

            rat_type det = abs(decomp_LU.determinant());

            return std::move(det);
        }

        template <typename rationalT>
        SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;
            using rat_vec = SparseVector<rationalT>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            return std::move(system_solve(b, decomp_LU));
        }

        template <typename rationalT, typename solverT>
        SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, solverT &solver)
        {
            using rat_type = rationalT;
            using rat_vec = SparseVector<rationalT>;
            using solver_type = solverT;

            return std::move(solver.solve(b));
        }
    } // namespace Utility

} // namespace ILPCuts
