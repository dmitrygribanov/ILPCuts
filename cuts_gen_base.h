#pragma once

#include "base_algo.h"

namespace ILPCuts
{
    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class cuts_gen_base : public dual_type_base_alg<intT, rationalT, lintT, lrationalT>
    {
    private:
        using my_type = cuts_gen_base<intT, rationalT, lintT, lrationalT>;
        using base = dual_type_base_alg<intT, rationalT, lintT, lrationalT>;

    public:
        DUAL_TYPE_ALGO_DEFS
    };
} // namespace ILPCuts
