#pragma once

#include "basic_utils.h"
#include <flint/fmpz_mat.h>

namespace ILPCuts
{
    namespace ModuleFlint
    {
        namespace Internal
        {
            template <typename int_type>
            void fill_mat(const SparseMatrix<int_type> &A, fmpz_mat_t A_flint);

            template <>
            void fill_mat(const SparseMatrix<int64_t> &A, fmpz_mat_t A_flint)
            {
                for (size_t k = 0; k < A.outerSize(); ++k)
                    for (typename SparseMatrix<int64_t>::InnerIterator it(A, k); it; ++it)
                    {
                        *fmpz_mat_entry(A_flint, it.row(), it.col()) = it.value();
                    }
            }

            template <>
            void fill_mat(const SparseMatrix<mpz_class> &A, fmpz_mat_t A_flint)
            {
                for (size_t k = 0; k < A.outerSize(); ++k)
                    for (typename SparseMatrix<mpz_class>::InnerIterator it(A, k); it; ++it)
                    {
                        fmpz *fmpz_entry = fmpz_mat_entry(A_flint, it.row(), it.col());
                        fmpz_set_mpz(fmpz_entry, it.value().get_mpz_t());
                    }
            }

            template <typename int_type>
            SparseMatrix<int_type> flint_to_Eigen(size_t m, size_t n, fmpz_mat_t A_flint)
            {
                using Triplet = Eigen::Triplet<int_type>;

                SparseMatrix<int_type> res(m, n);

                std::vector<Triplet> triplets_res;

                for (size_t i = 0; i < m; ++i)
                    for (size_t j = 0; j < n; ++j)
                    {
                        fmpz *flint_entry = fmpz_mat_entry(A_flint, i, j);
                        if (*flint_entry != 0)
                            triplets_res.push_back(Triplet(i, j, *flint_entry));
                    }

                res.setFromTriplets(triplets_res.begin(), triplets_res.end());
                res.makeCompressed();

                return res;
            }

            template <typename int_type>
            SparseMatrix<int_type> flint_to_Eigen(size_t n, fmpz_mat_t A_flint)
            {
                return flint_to_Eigen<int_type>(n, n, A_flint);
            }

            void
            fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A)
            {
                slong i, j, m, n;
                fmpz_mat_t A2, H2;
                flint_rand_t state;

                m = fmpz_mat_nrows(A);
                n = fmpz_mat_ncols(A);

                fmpz_mat_init(A2, m, n + m);
                fmpz_mat_init(H2, m, n + m);
                flint_randinit(state);

                /* append an identity matrix to the right of A */
                for (i = 0; i < m; i++)
                {
                    for (j = 0; j < n; j++)
                        fmpz_set(fmpz_mat_entry(A2, i, j), fmpz_mat_entry(A, i, j));
                    fmpz_one(fmpz_mat_entry(A2, i, n + i));
                }

                // fmpz_mat_hnf(H2, A2);
                // fmpz_mat_hnf_pernet_stein(H2, A2, state);
                fmpz_mat_hnf_xgcd(H2, A2);
                // fmpz_mat_hnf_classical(H2, A2);

                /* recover H and U */
                for (i = 0; i < m; i++)
                {
                    for (j = 0; j < n; j++)
                        fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H2, i, j));
                    for (j = n; j < n + m; j++)
                        fmpz_set(fmpz_mat_entry(U, i, j - n), fmpz_mat_entry(H2, i, j));
                }

                fmpz_mat_clear(A2);
                fmpz_mat_clear(H2);
                flint_randclear(state);
            }
        } // namespace Internal

        template <typename intT>
        void compute_HNF(const SparseMatrix<intT> &A, SparseMatrix<intT> &P, SparseMatrix<intT> &H);

        template <typename intT>
        void compute_HNF(const SparseMatrix<intT> &A, const SparseVector<intT> &x, SparseMatrix<intT> &H, SparseVector<intT> &y);

        template <typename intT>
        void compute_HNF(const SparseMatrix<intT> &A, SparseMatrix<intT> &H);

        template <typename intT>
        void compute_HNF(const SparseMatrix<intT> &A, const intT &det_mul, SparseMatrix<intT> &H);

        template <>
        void compute_HNF(const SparseMatrix<int64_t> &A, SparseMatrix<int64_t> &P, SparseMatrix<int64_t> &H)
        {
            size_t n = A.cols();
            size_t m = A.rows();

            fmpz_mat_t A_flint, H_flint, P_flint;
            fmpz_mat_init(A_flint, m, n);
            fmpz_mat_init(H_flint, m, n);
            fmpz_mat_init(P_flint, m, m);

            fmpz_mat_zero(A_flint);
            fmpz_mat_zero(H_flint);
            fmpz_mat_zero(P_flint);

            Internal::fill_mat(A, A_flint);

            // ILPCuts::Internal::fmpz_mat_hnf_transform(H_flint, P_flint, A_flint);
            fmpz_mat_hnf_transform(H_flint, P_flint, A_flint);

            P = Internal::flint_to_Eigen<int64_t>(m, P_flint);
            H = Internal::flint_to_Eigen<int64_t>(m, n, H_flint);

            fmpz_mat_clear(A_flint);
            fmpz_mat_clear(H_flint);
            fmpz_mat_clear(P_flint);
        }

        template <>
        void compute_HNF(const SparseMatrix<int64_t> &A, SparseMatrix<int64_t> &H)
        {
            size_t m = A.rows();
            size_t n = A.cols();

            fmpz_mat_t A_flint, H_flint;
            flint_rand_t state;
            fmpz_mat_init(A_flint, m, n);
            fmpz_mat_init(H_flint, m, n);
            flint_randinit(state);

            fmpz_mat_zero(A_flint);
            fmpz_mat_zero(H_flint);

            Internal::fill_mat(A, A_flint);

            fmpz_mat_hnf(H_flint, A_flint);
            // fmpz_mat_hnf_pernet_stein(H_flint, A_flint, state);

            H = Internal::flint_to_Eigen<int64_t>(m, n, H_flint);

            fmpz_mat_clear(A_flint);
            fmpz_mat_clear(H_flint);
            flint_randclear(state);
        }

        template <>
        void compute_HNF(const SparseMatrix<int64_t> &A, const SparseVector<int64_t> &x, SparseMatrix<int64_t> &H, SparseVector<int64_t> &y)
        {
            using int_vec = SparseMatrix<int64_t>;

            size_t m = A.rows();
            size_t n = A.cols();

            assert(m == x.size());

            int_vec A_ext = A;
            A_ext.conservativeResize(m, n + 1);

            // for (typename int_vec::InnerIterator it(x); it; ++it)
            // {
            //     A_ext.insert(it.index(), n) = it.value();
            // }

            A_ext.col(n) = x;
            A_ext.makeCompressed();

            // std::vector<size_t> cols = Utility::get_indices_of_nonzero_columns<int_type>(A_ext);

            // cols = Utility::get_indices_of_nonzero_columns<int_type>(A_ext);

            int_vec H_ext;
            compute_HNF(A_ext, H_ext);

            // cols = Utility::get_indices_of_nonzero_columns<int_type>(H_ext);

            H = H_ext.topLeftCorner(m, n);
            y = H_ext.col(n);
        }

        template <>
        void compute_HNF(const SparseMatrix<mpz_class> &A, const mpz_class &det_mul, SparseMatrix<mpz_class> &H)
        {
            size_t n = A.cols();
            fmpz_mat_t A_flint, H_flint;
            fmpz_t det_mul_flint;
            fmpz_mat_init(A_flint, n, n);
            fmpz_mat_init(H_flint, n, n);
            fmpz_init(det_mul_flint);

            fmpz_mat_zero(A_flint);
            fmpz_mat_zero(H_flint);

            Internal::fill_mat(A, A_flint);
            fmpz_set_mpz(det_mul_flint, det_mul.get_mpz_t());

            fmpz_mat_hnf_modular(H_flint, A_flint, det_mul_flint);

            H = Internal::flint_to_Eigen<mpz_class>(n, H_flint);

            fmpz_mat_clear(A_flint);
            fmpz_mat_clear(H_flint);
            fmpz_clear(det_mul_flint);
        }

        template <>
        void compute_HNF(const SparseMatrix<int64_t> &A, const int64_t &det_mul, SparseMatrix<int64_t> &H)
        {
            size_t n = A.cols();
            fmpz_mat_t A_flint, H_flint;
            fmpz_t det_mul_flint;
            fmpz_mat_init(A_flint, n, n);
            fmpz_mat_init(H_flint, n, n);
            fmpz_init(det_mul_flint);

            fmpz_mat_zero(A_flint);
            fmpz_mat_zero(H_flint);

            Internal::fill_mat(A, A_flint);
            fmpz_set_si(det_mul_flint, det_mul);

            fmpz_mat_hnf_modular(H_flint, A_flint, det_mul_flint);

            H = Internal::flint_to_Eigen<int64_t>(n, H_flint);

            fmpz_mat_clear(A_flint);
            fmpz_mat_clear(H_flint);
            fmpz_clear(det_mul_flint);
        }

        template <typename int_type>
        void compute_inverse(const SparseMatrix<int64_t> &A, int64_t &denum, SparseMatrix<int64_t> &Ai)
        {
            size_t n = A.cols();
            fmpz_mat_t A_flint, Ai_flint;
            fmpz_t denum_flint;
            fmpz_mat_init(A_flint, n, n);
            fmpz_mat_init(Ai_flint, n, n);
            fmpz_init(denum_flint);

            fmpz_mat_zero(A_flint);
            fmpz_mat_zero(Ai_flint);

            Internal::fill_mat(A, A_flint);

            fmpz_mat_inv(Ai_flint, denum_flint, A_flint);

            Ai = Internal::flint_to_Eigen<int64_t>(n, Ai_flint);
            denum = fmpz_get_si(denum_flint);

            fmpz_mat_clear(A_flint);
            fmpz_mat_clear(Ai_flint);
            fmpz_clear(denum_flint);
        }

        template <typename int_type>
        void compute_det(const SparseMatrix<int_type> &A, mpz_class &det);

        template <>
        void compute_det(const SparseMatrix<mpz_class> &A, mpz_class &det)
        {
            size_t n = A.cols();
            fmpz_mat_t A_flint;
            fmpz_t det_flint;
            fmpz_mat_init(A_flint, n, n);
            fmpz_init(det_flint);

            fmpz_mat_zero(A_flint);

            Internal::fill_mat<mpz_class>(A, A_flint);

            fmpz_mat_det(det_flint, A_flint);

            fmpz_get_mpz(det.get_mpz_t(), det_flint);

            fmpz_mat_clear(A_flint);
            fmpz_clear(det_flint);
        }
    } // namespace ModuleFlint
} // namespace ILPCuts
