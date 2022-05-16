#pragma once

#include "base_module.h"

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/det.h>
#include <linbox/solutions/solve.h>
#include <linbox/vector/sparse.h>

namespace ILPCuts
{
    namespace Utility
    {
        namespace LinBox
        {
            namespace Internal
            {
                template <typename Field>
                ::LinBox::DenseVector<Field> compute_some_kernel_vec(const Field &F, const ::LinBox::SparseMatrix<Field> &A)
                {
                    using dv_LB = ::LinBox::DenseVector<Field>;

                    size_t n = A.coldim();
                    dv_LB zero_vec_LB(F, n, F.zero);

                    dv_LB x_LB(F, n);
                    ::LinBox::Method::Blackbox method_LB;
                    ::LinBox::solve(x_LB, A, zero_vec_LB, method_LB);

                    return x_LB;
                }

                template <typename T>
                bool is_zero(const ::LinBox::DenseVector<T> &vec)
                {
                    for (size_t i = 0; i < vec.size(); ++i)
                    {
                        if (vec[i] != 0)
                            return false;
                    }

                    return true;
                }

                template <typename LinBoxValT, typename EigenValT>
                ::LinBox::SparseMatrix<LinBoxValT> create_from_Eigen(const Eigen::SparseMatrix<EigenValT> &A, const LinBoxValT &F)
                {
                    using eigen_mat_type = Eigen::SparseMatrix<EigenValT>;

                    size_t n = A.cols();
                    ::LinBox::SparseMatrix<LinBoxValT> intoA(F, n, n);

                    for (size_t k = 0; k < A.outerSize(); ++k)
                        for (typename eigen_mat_type::InnerIterator it(A, k); it; ++it)
                        {
                            typename LinBoxValT::Element new_F_elem;
                            F.init(new_F_elem, it.value());
                            intoA.setEntry(it.row(), it.col(), new_F_elem);
                        }

                    return intoA;
                }

            }

            template <typename int_type>
            int_type compute_det(Eigen::SparseMatrix<int_type> A);

            template <>
            mpz_class compute_det(Eigen::SparseMatrix<mpz_class> A)
            {
                using eigen_mat_type = Eigen::SparseMatrix<mpz_class>;

                prune_zeroes_int(A);

                using ZRing = Givaro::ZRing<Givaro::Integer>;
                ZRing Z;

                // using QField = Givaro::QField<Givaro::Rational>;
                // QField Q;
                ::LinBox::SparseMatrix<ZRing> A_LinBox = Internal::create_from_Eigen<ZRing, mpz_class>(A, Z);

                ZRing::Element det_LinBox;
                ::LinBox::Method::Blackbox method_LinBox;
                ::LinBox::det(det_LinBox, A_LinBox, method_LinBox);

                return mpz_class(det_LinBox.get_mpz_const());
            }

            template <typename int_type>
            Eigen::SparseVector<int_type> compute_some_kernel_vec(Eigen::SparseMatrix<int_type> A, const int_type &p);

            template <>
            Eigen::SparseVector<int64_t> compute_some_kernel_vec(Eigen::SparseMatrix<int64_t> A, const int64_t &p)
            {
                using Field = Givaro::Modular<int64_t>;
                using eigen_mat_type = Eigen::SparseMatrix<int64_t>;

                A.unaryExpr(Utility::mod_unaryop(p));
                prune_zeroes_int(A);

                size_t n = A.cols();
                Field F(p);
                ::LinBox::SparseMatrix<Field> A_LB = Internal::create_from_Eigen<Field, int64_t>(A, F);

                ::LinBox::DenseVector<Field> res_LB;
                bool has_non_zero_solution = false;

                while (!has_non_zero_solution)
                {
                    res_LB = Internal::compute_some_kernel_vec(F, A_LB);

                    has_non_zero_solution = !Internal::is_zero<Field>(res_LB);
                }

                Eigen::SparseVector<int64_t> res(n);
                for (size_t i = 0; i < n; ++i)
                {
                    if (res_LB[i] != 0)
                        res.insert(i) = res_LB[i];
                }

                return res;
            }

            template <>
            Eigen::SparseVector<mpz_class> compute_some_kernel_vec(Eigen::SparseMatrix<mpz_class> A, const mpz_class &p)
            {
                using Field = Givaro::Modular<Givaro::Integer>;
                using eigen_mat_type = Eigen::SparseMatrix<mpz_class>;

                A.unaryExpr(Utility::mod_unaryop(p));
                prune_zeroes_int(A);
                size_t n = A.cols();
                Field F(p);
                ::LinBox::SparseMatrix<Field> A_LB = Internal::create_from_Eigen<Field, mpz_class>(A, F);

                ::LinBox::DenseVector<Field> res_LB;
                bool has_non_zero_solution = false;

                while (!has_non_zero_solution)
                {
                    res_LB = Internal::compute_some_kernel_vec(F, A_LB);
                    has_non_zero_solution = !Internal::is_zero<Field>(res_LB);
                }

                Eigen::SparseVector<mpz_class> res(n);
                for (size_t i = 0; i < n; ++i)
                {
                    if (res_LB[i] != 0)
                        res.insert(i) = mpz_class(res_LB[i].get_mpz_const());
                }

                return res;
            }

            template <typename int_type>
            Eigen::SparseVector<int_type> solve_ZZ(Eigen::SparseMatrix<int_type> A, const Eigen::SparseVector<int_type> &b);

            template <>
            Eigen::SparseVector<mpz_class> solve_ZZ(Eigen::SparseMatrix<mpz_class> A, const Eigen::SparseVector<mpz_class> &b)
            {
                prune_zeroes_int(A);
                size_t n = A.cols();

                using ZRing = Givaro::ZRing<Givaro::Integer>;
                ZRing Z;

                ::LinBox::SparseMatrix<ZRing> A_LB = Internal::create_from_Eigen<ZRing, mpz_class>(A, Z);
            }
        }
    }
}