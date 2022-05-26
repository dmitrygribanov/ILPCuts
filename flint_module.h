#pragma once

#include "base_module.h"
#include "flint/fmpz_mat.h"

namespace ILPCuts
{
    namespace Utility
    {
        namespace flint
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
                SparseMatrix<int_type> flint_to_Eigen(size_t n, fmpz_mat_t A_flint)
                {
                    SparseMatrix<int_type> res(n, n);

                    for (size_t i = 0; i < n; ++i)
                        for (size_t j = 0; j < n; ++j)
                        {
                            fmpz *flint_entry = fmpz_mat_entry(A_flint, i, j);
                            if (*flint_entry != 0)
                                res.insert(i, j) = *flint_entry;
                        }

                    return res;
                }
            } // namespace Internal

            template <typename int_type>
            void compute_HNF(const SparseMatrix<int_type> &A, SparseMatrix<int_type> &P, SparseMatrix<int_type> &H);

            template <typename int_type>
            void compute_HNF(const SparseMatrix<int_type> &A, SparseMatrix<int_type> &H);

            template <typename int_type>
            void compute_HNF(const SparseMatrix<int_type> &A, const int_type &det_mul, SparseMatrix<int_type> &H);

            template <>
            void compute_HNF(const SparseMatrix<int64_t> &A, SparseMatrix<int64_t> &P, SparseMatrix<int64_t> &H)
            {
                size_t n = A.cols();
                fmpz_mat_t A_flint, H_flint, P_flint;
                fmpz_mat_init(A_flint, n, n);
                fmpz_mat_init(H_flint, n, n);
                fmpz_mat_init(P_flint, n, n);

                fmpz_mat_zero(A_flint);
                fmpz_mat_zero(H_flint);
                fmpz_mat_zero(P_flint);

                Internal::fill_mat(A, A_flint);

                fmpz_mat_hnf_transform(H_flint, P_flint, A_flint);

                P = Internal::flint_to_Eigen<int64_t>(n, P_flint);
                H = Internal::flint_to_Eigen<int64_t>(n, H_flint);

                fmpz_mat_clear(A_flint);
                fmpz_mat_clear(H_flint);
                fmpz_mat_clear(P_flint);
            }

            template <>
            void compute_HNF(const SparseMatrix<int64_t> &A, SparseMatrix<int64_t> &H)
            {
                size_t n = A.cols();
                fmpz_mat_t A_flint, H_flint;
                fmpz_mat_init(A_flint, n, n);
                fmpz_mat_init(H_flint, n, n);

                fmpz_mat_zero(A_flint);
                fmpz_mat_zero(H_flint);

                Internal::fill_mat(A, A_flint);

                fmpz_mat_hnf(H_flint, A_flint);

                H = Internal::flint_to_Eigen<int64_t>(n, H_flint);

                fmpz_mat_clear(A_flint);
                fmpz_mat_clear(H_flint);
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
        } // namespace flint
    }     // namespace Utility

} // namespace ILPCuts
