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

        // template <>
        // void prune_zeroes_int<zz_p>(SparseVector<zz_p> &A)
        // {
        //     A.prune(KeepIntFunc<zz_p>());
        // }
    } // namespace Utility
}