#pragma once

#include "conic_cuts_base.h"

namespace ILPCuts
{
    template <typename intT, typename rationalT, typename lintT, typename lratT>
    class ConicCutsGenerator_Zp : public ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>
    {
    private:
        using base = ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>;

    public:
        using typename base::int_type;
        using typename base::lint_type;
        using typename base::lrational_type;
        using typename base::rational_type;

        using typename base::is_mat;
        using typename base::lis_mat;
        using typename base::lrs_mat;
        using typename base::rs_mat;

        using typename base::is_vec;
        using typename base::lis_vec;
        using typename base::lrs_vec;
        using typename base::rs_vec;

        using typename base::cone_traits;
        using typename base::dir_coeffs;
        using typename base::lsparce_transform;
        using typename base::process_func;
        using typename base::sparce_transform;
        using typename base::transform_base;
    };
}