#pragma once

#include "basic_utils.h"

#include <sstream>
#include <NTL/mat_GF2.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/HNF.h>

namespace ILPCuts
{
    namespace Utility
    {
        namespace NTL
        {
            using namespace ::NTL;

            namespace Internal
            {
                ZZ mpz_to_ZZ(const mpz_class &a)
                {
                    std::string str = a.get_str();
                    return ::NTL::conv<ZZ>(str.c_str());
                }

                mpz_class ZZ_to_mpz(const ZZ &a)
                {
                    std::stringstream ss;
                    ss << a;
                    return mpz_class(ss.str());
                }

                template <typename in_type, typename out_type>
                struct Conversion
                {
                    out_type operator()(const in_type &a) const { return out_type(a); }
                };

                template <>
                struct Conversion<ZZ, int64_t>
                {
                    int64_t operator()(const ZZ &a) const { return conv<int64_t>(a); }
                };

                template <>
                struct Conversion<mpz_class, ZZ_p>
                {
                    ZZ_p operator()(const mpz_class &a) const { return conv<ZZ_p>(mpz_to_ZZ(a)); }
                };

                template <>
                struct Conversion<ZZ_p, mpz_class>
                {
                    mpz_class operator()(const ZZ_p &a) const { return ZZ_to_mpz(conv<ZZ>(a)); }
                };

                template <>
                struct Conversion<mpz_class, GF2>
                {
                    GF2 operator()(const mpz_class &a) const { return ::NTL::conv<GF2>(a.get_si()); }
                };

                template <>
                struct Conversion<GF2, mpz_class>
                {
                    mpz_class operator()(const GF2 &a) const { return mpz_class(conv<int>(a)); }
                };

                template <>
                struct Conversion<zz_p, int64_t>
                {
                    int64_t operator()(const zz_p &a) const { return conv<int64_t>(a); }
                };

                template <>
                struct Conversion<GF2, int64_t>
                {
                    int64_t operator()(const GF2 &a) const { return conv<int64_t>(a); }
                };

                template <>
                struct Conversion<ZZ, mpz_class>
                {
                    mpz_class operator()(const ZZ &a) const { return ZZ_to_mpz(a); }
                };

                template <>
                struct Conversion<mpz_class, ZZ>
                {
                    ZZ operator()(const mpz_class &a) const { return mpz_to_ZZ(a); }
                };

                template <typename in_type, typename out_type, typename conversion_type = Conversion<in_type, out_type>>
                Mat<out_type> eigen_to_ntl(const Eigen::SparseMatrix<in_type> &A, conversion_type conv = conversion_type())
                {
                    size_t n = A.cols();
                    Mat<out_type> res;
                    res.SetDims(n, n);

                    for (size_t k = 0; k < A.outerSize(); ++k)
                        for (typename Eigen::SparseMatrix<in_type>::InnerIterator it(A, k); it; ++it)
                        {
                            if (it.value() != 0)
                                res[it.row()][it.col()] = conv(it.value());
                        }

                    return res;
                }

                template <typename in_type, typename out_type, typename conversion_type = Conversion<in_type, out_type>>
                Vec<out_type> eigen_to_ntl(const Eigen::SparseVector<in_type> &v, conversion_type conv = conversion_type())
                {
                    size_t n = v.size();
                    Vec<out_type> res;
                    res.SetLength(n);

                    for (typename Eigen::SparseVector<in_type>::InnerIterator it(v); it; ++it)
                    {
                        if (it.value() != 0)
                            res(it.index()) = conv(it.value());
                    }

                    return res;
                }

                template <typename in_type, typename out_type, typename conversion_type = Conversion<in_type, out_type>>
                Eigen::SparseMatrix<out_type> ntl_to_eigen(const Mat<in_type> &A, conversion_type conv = conversion_type())
                {
                    size_t n = A.NumCols();
                    Eigen::SparseMatrix<out_type> res(n, n);

                    for (size_t i = 0; i < n; ++i)
                        for (size_t j = 0; j < n; ++j)
                            if (A[i][j] != 0)
                                res.insert(i, j) = conv(A[i][j]);

                    return res;
                }

                template <typename in_type, typename out_type, typename conversion_type = Conversion<in_type, out_type>>
                Eigen::SparseVector<out_type> ntl_to_eigen(const Vec<in_type> &v, conversion_type conv = conversion_type())
                {
                    size_t n = v.length();
                    Eigen::SparseVector<out_type> res(n);

                    for (size_t i = 0; i < n; ++i)
                        if (v[i] != 0)
                            res.insert(i) = conv(v[i]);

                    return res;
                }

                template <typename T>
                size_t sparsest_row_index(const Mat<T> &A)
                {
                    size_t n = A.NumCols();
                    size_t m = A.NumRows();
                    size_t i_min_nonzeros = 0;
                    size_t min_nonzeros = n;

                    for (size_t i = 0; i < m; ++i)
                    {
                        size_t nonzeros = 0;
                        for (size_t j = 0; j < n; ++j)
                        {
                            if (A[i][j] != 0)
                                ++nonzeros;
                        }
                        if (nonzeros < min_nonzeros)
                        {
                            min_nonzeros = nonzeros;
                            i_min_nonzeros = i;
                        }
                    }

                    return i_min_nonzeros;
                }

            }

            template <typename int_type>
            Eigen::SparseVector<int_type> kernel_vec_GF2(Eigen::SparseMatrix<int_type> A)
            {
                size_t n = A.cols();

                A = A.unaryExpr(std::ref(Utility::mod2_unaryop<int_type>)).eval();
                prune_zeroes_int(A);
                A = A.transpose().eval();

                mat_GF2 A_ntl = Internal::eigen_to_ntl<int_type, GF2>(A);
                mat_GF2 X_ntl;

                kernel(X_ntl, A_ntl);

                size_t sprsest_row_idx = Internal::sparsest_row_index(X_ntl);

                return Internal::ntl_to_eigen<GF2, int_type>(X_ntl[sprsest_row_idx]);
            }

            template <typename int_type>
            Eigen::SparseVector<int_type> kernel_vec_Zp(Eigen::SparseMatrix<int_type> A, const int_type &p);

            template <>
            Eigen::SparseVector<int64_t> kernel_vec_Zp(Eigen::SparseMatrix<int64_t> A, const int64_t &p)
            {
                size_t n = A.cols();

                zz_p::init(p);

                A = A.unaryExpr(mod_unaryop(p)).eval();
                prune_zeroes_int(A);
                A = A.transpose().eval();

                mat_zz_p A_ntl = Internal::eigen_to_ntl<int64_t, zz_p>(A);
                mat_zz_p X_ntl;

                kernel(X_ntl, A_ntl);

                size_t sprsest_row_idx = Internal::sparsest_row_index(X_ntl);

                return Internal::ntl_to_eigen<zz_p, int64_t>(X_ntl[sprsest_row_idx]);
            }

            template <>
            Eigen::SparseVector<mpz_class> kernel_vec_Zp(Eigen::SparseMatrix<mpz_class> A, const mpz_class &p)
            {
                size_t n = A.cols();

                ZZ_p::init(Internal::mpz_to_ZZ(p));

                A = A.unaryExpr(mod_unaryop(p)).eval();
                prune_zeroes_int(A);
                A = A.transpose().eval();

                mat_ZZ_p A_ntl = Internal::eigen_to_ntl<mpz_class, ZZ_p>(A);
                mat_ZZ_p X_ntl;

                kernel(X_ntl, A_ntl);

                size_t sprsest_row_idx = Internal::sparsest_row_index(X_ntl);

                return Internal::ntl_to_eigen<ZZ_p, mpz_class>(X_ntl[sprsest_row_idx]);
            }

            // template <>
            // Eigen::SparseVector<mpz_class> kernel_vec_Zp(Eigen::SparseMatrix<mpz_class> A, const mpz_class &p)
            // {
            //     size_t n = A.cols();

            //     A.unaryExpr(mod_unaryop(p));
            //     prune_zeroes_int(A);

            //     mat_ZZ_p A_ntl = Internal::eigen_to_ntl<mpz_class, ZZ_p>(A);
            //     mat_ZZ_p X_ntl;

            //     kernel(X_ntl, A_ntl);

            //     return Internal::ntl_to_eigen<ZZ_p, mpz_class>(X_ntl);
            // }

            template <typename int_type>
            Eigen::SparseVector<int_type> solve_ZZ(Eigen::SparseMatrix<int_type> A, Eigen::SparseVector<int_type> b)
            {
                size_t n = A.cols();

                prune_zeroes_int(A);
                A = A.transpose().eval();

                vec_ZZ b_ntl = Internal::eigen_to_ntl<int_type, ZZ>(b);
                mat_ZZ A_ntl = Internal::eigen_to_ntl<int_type, ZZ>(A);

                vec_ZZ x_ntl;
                ZZ d;
                solve(d, x_ntl, A_ntl, b_ntl);

                return Internal::ntl_to_eigen<ZZ, int_type>(x_ntl);
            }

            template <typename int_type>
            int_type det(Eigen::SparseMatrix<int_type> A)
            {
                mat_ZZ A_ntl = Internal::eigen_to_ntl<int_type, ZZ>(A);

                ZZ d;
                determinant(d, A_ntl);

                return Internal::Conversion<ZZ, int_type>()(d);
            }

            template <typename int_type>
            std::tuple<SparseMatrix<int_type>, int_type> compute_inverse(const SparseMatrix<int_type> &A);

            template <>
            std::tuple<SparseMatrix<mpz_class>, mpz_class> compute_inverse(const SparseMatrix<mpz_class> &A)
            {
                mat_ZZ A_ntl = Internal::eigen_to_ntl<mpz_class, ZZ>(A);

                ZZ d_ntl;
                mat_ZZ Ai_ntl;
                inv(d_ntl, Ai_ntl, A_ntl);

                return std::make_tuple(Internal::ntl_to_eigen<ZZ, mpz_class>(Ai_ntl), Internal::ZZ_to_mpz(d_ntl));
            }

            template <typename int_type>
            SparseMatrix<int_type> compute_HNF(const SparseMatrix<int_type> &A, const int_type &det_mul)
            {
                mat_ZZ A_ntl = Internal::eigen_to_ntl<int_type, ZZ>(A);

                ZZ det_mul_ntl = Internal::Conversion<int_type, ZZ>()(det_mul);

                mat_ZZ H_ntl;
                HNF(H_ntl, A_ntl, det_mul_ntl);

                return Internal::ntl_to_eigen<ZZ, int_type>(H_ntl);
            }
        }
    }
}