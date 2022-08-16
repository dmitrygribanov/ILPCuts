#pragma once

#include <NTL/lzz_p.h>

#include <eigen3/Eigen/Core>

#include "basic_utils.h"

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

    namespace Utility
    {
        template <>
        void prune_zeroes_int<zz_p>(SparseMatrix<zz_p> &A)
        {
            A.prune(KeepIntFunc<zz_p>());
        }

        template <>
        struct conv<zz_p, int64_t>
        {
            static int64_t make(const zz_p &a)
            {
                return rep(a);
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

    } // namespace Utility
}