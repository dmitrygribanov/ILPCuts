#pragma once

#include "basic_utils.h"

#include <NTL/lzz_p.h>
#include <NTL/matrix.h>
#include <NTL/mat_GF2.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/HNF.h>

#include <eigen3/Eigen/Core>

namespace Eigen
{
    template <>
    struct NumTraits<NTL::zz_p> : GenericNumTraits<NTL::zz_p>
    {
        // typedef NTL::zz_p Real;
        // typedef NTL::zz_p NonInteger;
        // typedef NTL::zz_p Nested;

        // static inline Real epsilon() { return 0; }
        // static inline Real dummy_precision() { return 0; }
        // static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 1,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 50,
            MulCost = 75
        };
    };
}

namespace ILPCuts
{
    using namespace NTL;
    namespace ModuleNTL
    {
        ZZ mpz_to_ZZ(const mpz_class &a)
        {
            std::string str = a.get_str();
            return NTL::conv<ZZ>(str.c_str());
        }

        mpz_class ZZ_to_mpz(const ZZ &a)
        {
            std::stringstream ss;
            ss << a;
            return mpz_class(ss.str());
        }
    } // namespace ModuleNTL

    namespace Utility
    {
        template <>
        void prune_zeroes_int<zz_p>(SparseMatrix<zz_p> &A)
        {
            A.prune(KeepIntFunc<zz_p>());
        }
    } // namespace Utility

    template <>
    struct conv<zz_p, int64_t>
    {
        static int64_t make(const zz_p &a)
        {
            return rep(a);
        }
    };

    template <>
    struct conv<ZZ, int64_t>
    {
        static int64_t make(const ZZ &a)
        {
            return NTL::conv<int64_t>(a);
        }
    };

    template <>
    struct conv<mpz_class, ZZ_p>
    {
        static ZZ_p make(const mpz_class &a)
        {
            return NTL::conv<ZZ_p>(ModuleNTL::mpz_to_ZZ(a));
        }
    };

    template <>
    struct conv<ZZ_p, mpz_class>
    {
        static mpz_class make(const ZZ_p &a)
        {
            return ModuleNTL::ZZ_to_mpz(NTL::conv<ZZ>(a));
        }
    };

    template <>
    struct conv<mpz_class, GF2>
    {
        static GF2 make(const mpz_class &a)
        {
            return NTL::conv<GF2>(a.get_si());
        }
    };

    template <>
    struct conv<GF2, mpz_class>
    {
        static mpz_class make(const GF2 &a)
        {
            return mpz_class(NTL::conv<int>(a));
        }
    };

    template <>
    struct conv<GF2, int64_t>
    {
        static int64_t make(const GF2 &a)
        {
            return NTL::conv<int64_t>(a);
        }
    };

    template <>
    struct conv<ZZ, mpz_class>
    {
        static mpz_class make(const ZZ &a)
        {
            return ModuleNTL::ZZ_to_mpz(a);
        }
    };

    template <>
    struct conv<mpz_class, ZZ>
    {
        static ZZ make(const mpz_class &a)
        {
            return ModuleNTL::mpz_to_ZZ(a);
        }
    };

    template <typename inT, typename outT>
    struct conv<NTL::Vec<inT>, SparseVector<outT>>
    {
        static SparseVector<outT> make(const NTL::Vec<inT> &a)
        {
            size_t n = a.length();

            SparseVector<outT> res(n);
            for (size_t i = 0; i < n; ++i)
            {
                if (a[i] != 0)
                    res.insertBack(i) = conv<inT, outT>::make(a[i]);
            }

            return res;
        }
    };

    template <typename inT, typename outT>
    struct conv<SparseVector<inT>, NTL::Vec<outT>>
    {
        static NTL::Vec<outT> make(const SparseVector<inT> &a)
        {
            size_t n = a.size();

            NTL::Vec<outT> res;
            res.SetLength(n, outT(0));

            for (typename SparseVector<inT>::InnerIterator it(a); it; ++it)
            {
                res[it.index()] = conv<inT, outT>::make(it.value());
            }

            return res;
        }
    };

    template <typename inT, typename outT>
    struct conv<NTL::Mat<inT>, SparseMatrix<outT>>
    {
        static SparseMatrix<outT> make(const NTL::Mat<inT> &A)
        {
            size_t n = A.NumCols();
            size_t m = A.NumRows();
            Eigen::SparseMatrix<outT> res(m, n);

            for (size_t i = 0; i < m; ++i)
                for (size_t j = 0; j < n; ++j)
                    if (A[i][j] != 0)
                        res.insert(i, j) = conv<inT, outT>::make(A[i][j]);

            return res;
        }
    };

    template <typename inT, typename outT>
    struct conv<SparseMatrix<inT>, NTL::Mat<outT>>
    {
        static NTL::Mat<outT> make(const SparseMatrix<inT> &A)
        {
            size_t n = A.cols();
            size_t m = A.rows();
            Mat<outT> res;
            res.SetDims(m, n);

            for (size_t k = 0; k < A.outerSize(); ++k)
                for (typename Eigen::SparseMatrix<inT>::InnerIterator it(A, k); it; ++it)
                {
                    if (it.value() != 0)
                        res[it.row()][it.col()] = conv<inT, outT>::make(it.value());
                }

            return res;
        }
    };
    namespace ModuleNTL
    {
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

        template <typename int_type>
        Eigen::SparseVector<int_type> kernel_vec_GF2(Eigen::SparseMatrix<int_type> A)
        {
            size_t n = A.cols();

            A = A.unaryExpr(std::ref(Utility::mod2_unaryop<int_type>)).eval();
            Utility::prune_zeroes_int(A);
            A = A.transpose().eval();

            mat_GF2 A_ntl = conv<Eigen::SparseMatrix<int_type>, mat_GF2>::make(A);
            mat_GF2 X_ntl;

            NTL::kernel(X_ntl, A_ntl);

            size_t sprsest_row_idx = sparsest_row_index(X_ntl);

            return conv<vec_GF2, Eigen::SparseVector<int_type>>::make(X_ntl[sprsest_row_idx]);
        }

        template <typename int_type>
        Eigen::SparseVector<int_type> kernel_vec_Zp(Eigen::SparseMatrix<int_type> A, const int_type &p);

        template <>
        Eigen::SparseVector<int64_t> kernel_vec_Zp(Eigen::SparseMatrix<int64_t> A, const int64_t &p)
        {
            size_t n = A.cols();

            zz_p::init(p);

            A = A.unaryExpr(Utility::mod_unaryop(p)).eval();
            Utility::prune_zeroes_int(A);
            A = A.transpose().eval();

            mat_zz_p A_ntl = conv<Eigen::SparseMatrix<int64_t>, mat_zz_p>::make(A);
            mat_zz_p X_ntl;

            NTL::kernel(X_ntl, A_ntl);

            size_t sprsest_row_idx = sparsest_row_index(X_ntl);

            return conv<vec_zz_p, Eigen::SparseVector<int64_t>>::make(X_ntl[sprsest_row_idx]);
        }

        template <>
        Eigen::SparseVector<mpz_class> kernel_vec_Zp(Eigen::SparseMatrix<mpz_class> A, const mpz_class &p)
        {
            size_t n = A.cols();

            ZZ_p::init(conv<mpz_class, ZZ>::make(p));

            A = A.unaryExpr(Utility::mod_unaryop(p)).eval();
            Utility::prune_zeroes_int(A);
            A = A.transpose().eval();

            mat_ZZ_p A_ntl = conv<Eigen::SparseMatrix<mpz_class>, mat_ZZ_p>::make(A);
            mat_ZZ_p X_ntl;

            kernel(X_ntl, A_ntl);

            size_t sprsest_row_idx = sparsest_row_index(X_ntl);

            return conv<vec_ZZ_p, Eigen::SparseVector<mpz_class>>::make(X_ntl[sprsest_row_idx]);
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

            Utility::prune_zeroes_int(A);
            A = A.transpose().eval();

            vec_ZZ b_ntl = conv<Eigen::SparseVector<int_type>, vec_ZZ>::make(b);
            mat_ZZ A_ntl = conv<Eigen::SparseMatrix<int_type>, mat_ZZ>::make(A);

            vec_ZZ x_ntl;
            ZZ d;
            NTL::solve(d, x_ntl, A_ntl, b_ntl);

            return conv<vec_ZZ, Eigen::SparseVector<int_type>>::make(x_ntl);
        }

        template <typename int_type>
        int_type det(Eigen::SparseMatrix<int_type> A)
        {
            mat_ZZ A_ntl = conv<Eigen::SparseMatrix<int_type>, mat_ZZ>::make(A);

            ZZ d;
            NTL::determinant(d, A_ntl);

            return conv<ZZ, int_type>::make(d);
        }

        template <typename int_type>
        std::tuple<SparseMatrix<int_type>, int_type> compute_inverse(const SparseMatrix<int_type> &A);

        template <>
        std::tuple<SparseMatrix<mpz_class>, mpz_class> compute_inverse(const SparseMatrix<mpz_class> &A)
        {
            mat_ZZ A_ntl = conv<SparseMatrix<mpz_class>, mat_ZZ>::make(A);

            ZZ d_ntl;
            mat_ZZ Ai_ntl;
            NTL::inv(d_ntl, Ai_ntl, A_ntl);

            return make_tuple(conv<mat_ZZ, SparseMatrix<mpz_class>>::make(Ai_ntl), conv<ZZ, mpz_class>::make(d_ntl));
        }

        template <typename int_type>
        SparseMatrix<int_type> compute_HNF(const SparseMatrix<int_type> &A, const int_type &det_mul)
        {
            mat_ZZ A_ntl = conv<SparseMatrix<int_type>, mat_ZZ>::make(A);

            ZZ det_mul_ntl = conv<int_type, ZZ>::make(det_mul);

            mat_ZZ H_ntl;
            NTL::HNF(H_ntl, A_ntl, det_mul_ntl);

            return conv<mat_ZZ, SparseMatrix<int_type>>::make(H_ntl);
        }

    } // namespace ModuleNTL
} // namespace ILPCuts