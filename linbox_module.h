#pragma once

#include "basic_utils.h"

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
                std::tuple<::LinBox::DenseVector<Field>, Field> solve(const Field &F, const ::LinBox::SparseMatrix<Field> &A, const ::LinBox::DenseVector<Field> &b)
                {
                    using dv_LB_type = ::LinBox::DenseVector<Field>;

                    size_t n = A.coldim();

                    dv_LB_type x_LB(F, n);
                    typename Field::Element denum;
                    ::LinBox::Method::Blackbox method_LB;
                    method_LB.preconditioner = ::LinBox::Preconditioner::Sparse;
                    method_LB.singularity = ::LinBox::Singularity::NonSingular;
                    ::LinBox::solve(x_LB, denum, A, b, method_LB);

                    return std::make_tuple(x_LB, denum);
                }

                template <typename T>
                bool is_zero(const ::LinBox::DenseVector<T> &v)
                {
                    for (size_t i = 0; i < v.size(); ++i)
                    {
                        if (v[i] != 0)
                            return false;
                    }

                    return true;
                }

                template <typename Field>
                ::LinBox::DenseVector<Field> compute_some_kernel_vec(const Field &F, const ::LinBox::SparseMatrix<Field> &A)
                {
                    using dv_LB_type = ::LinBox::DenseVector<Field>;

                    size_t n = A.coldim();
                    dv_LB_type zero_vec_LB(F, n, F.zero);

                    ::LinBox::DenseVector<Field> x_LB(F, n);
                    bool has_non_zero_solution = false;
                    while (!has_non_zero_solution)
                    {
                        ::LinBox::Method::SparseElimination method_LB;
                        method_LB.preconditioner = ::LinBox::Preconditioner::Sparse;
                        method_LB.singularity = ::LinBox::Singularity::Singular;
                        ::LinBox::solve(x_LB, A, zero_vec_LB, method_LB);
                        has_non_zero_solution = !Internal::is_zero<Field>(x_LB);
                    }

                    return x_LB;
                }

                template <typename in_type, typename out_type>
                struct Conversion
                {
                    out_type operator()(const typename in_type::Element &a) const { return out_type(a); }
                };

                template <>
                struct Conversion<Givaro::Modular<Givaro::Integer>, mpz_class>
                {
                    mpz_class operator()(const typename Givaro::Modular<Givaro::Integer>::Element &a) const { return mpz_class(a.get_mpz_const()); }
                };
                template <>
                struct Conversion<Givaro::ZRing<Givaro::Integer>, mpz_class>
                {
                    mpz_class operator()(const typename Givaro::ZRing<Givaro::Integer>::Element &a) const { return mpz_class(a.get_mpz_const()); }
                };

                template <typename inT, typename outT, typename ConversionT = Conversion<inT, outT>>
                Eigen::SparseVector<outT> Eigen_from_LinBox(const ::LinBox::DenseVector<inT> &v, ConversionT conv = ConversionT())
                {
                    size_t n = v.size();
                    Eigen::SparseVector<outT> res(n);
                    for (size_t i = 0; i < n; ++i)
                    {
                        if (v[i] != 0)
                            res.insert(i) = conv(v[i]);
                    }

                    return res;
                }

                template <typename intT, typename outT>
                ::LinBox::SparseMatrix<outT> LinBox_from_Eigen(const Eigen::SparseMatrix<intT> &A, const outT &F)
                {
                    using eigen_mat_type = Eigen::SparseMatrix<intT>;

                    size_t n = A.cols();
                    ::LinBox::SparseMatrix<outT> res_A(F, n, n);

                    for (size_t k = 0; k < A.outerSize(); ++k)
                        for (typename eigen_mat_type::InnerIterator it(A, k); it; ++it)
                        {
                            typename outT::Element new_F_elem;
                            F.init(new_F_elem, it.value());
                            res_A.setEntry(it.row(), it.col(), new_F_elem);
                        }

                    return res_A;
                }

                template <typename intT, typename outT>
                ::LinBox::DenseVector<outT> LinBox_from_Eigen(const Eigen::SparseVector<intT> &vec, const outT &F)
                {
                    using eigen_vec_type = Eigen::SparseVector<intT>;

                    size_t n = vec.size();
                    ::LinBox::DenseVector<outT> res_vec(F, n);

                    for (typename eigen_vec_type::InnerIterator it(vec); it; ++it)
                    {
                        typename outT::Element new_F_elem;
                        F.init(new_F_elem, it.value());
                        res_vec.setEntry(it.index(), new_F_elem);
                    }

                    return res_vec;
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
                ::LinBox::SparseMatrix<ZRing> A_LinBox = Internal::LinBox_from_Eigen<mpz_class, ZRing>(A, Z);

                ZRing::Element det_LinBox;
                ::LinBox::Method::Blackbox method_LinBox;
                method_LinBox.singularity = ::LinBox::Singularity::NonSingular;
                ::LinBox::det(det_LinBox, A_LinBox, method_LinBox);

                return mpz_class(det_LinBox.get_mpz_const());
            }

            template <typename int_type>
            Eigen::SparseVector<int_type> compute_some_kernel_vec(Eigen::SparseMatrix<int_type> A, const int_type &p);

            template <>
            Eigen::SparseVector<int64_t> compute_some_kernel_vec(Eigen::SparseMatrix<int64_t> A, const int64_t &p)
            {
                using Field = Givaro::Modular<int64_t>;

                A = A.unaryExpr(Utility::mod_unaryop(p)).eval();
                prune_zeroes_int(A);

                size_t n = A.cols();
                Field F(p);
                ::LinBox::SparseMatrix<Field> A_LB = Internal::LinBox_from_Eigen<int64_t, Field>(A, F);

                ::LinBox::DenseVector<Field> res_LB = Internal::compute_some_kernel_vec(F, A_LB);

                return Internal::Eigen_from_LinBox<Field, int64_t>(res_LB);
            }

            template <>
            Eigen::SparseVector<mpz_class> compute_some_kernel_vec(Eigen::SparseMatrix<mpz_class> A, const mpz_class &p)
            {
                using Field = Givaro::Modular<Givaro::Integer>;

                A = A.unaryExpr(Utility::mod_unaryop(p)).eval();
                prune_zeroes_int(A);

                size_t n = A.cols();
                Field F(p);
                ::LinBox::SparseMatrix<Field> A_LB = Internal::LinBox_from_Eigen<mpz_class, Field>(A, F);

                ::LinBox::DenseVector<Field> res_LB = Internal::compute_some_kernel_vec(F, A_LB);

                return Internal::Eigen_from_LinBox<Field, mpz_class>(res_LB);
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

                ::LinBox::SparseMatrix<ZRing> A_LB = Internal::LinBox_from_Eigen<mpz_class, ZRing>(A, Z);
                ::LinBox::DenseVector<ZRing> b_LB = Internal::LinBox_from_Eigen<mpz_class, ZRing>(b, Z);

                auto [x_LB, denum] = Internal::solve(Z, A_LB, b_LB);

                assert(denum == 1);

                return Internal::Eigen_from_LinBox<ZRing, mpz_class>(x_LB);
            }
        }
    }
}