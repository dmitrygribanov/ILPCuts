#pragma once

#include <numeric>
#include <functional>
#include <memory>

#include "basic_utils.h"
#include "gmp_supp.h"
#include "timer.h"
#include "cuts_strategies.h"
#include "sparse_transform.h"
#include "boost_tee_logging.h"

namespace ILPCuts
{
    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class base_alg
    {
    public:
        using rational_type = rationalT;
        using int_type = intT;

        using lrational_type = lrationalT;
        using lint_type = lintT;

        using is_mat = SparseMatrix<int_type>;
        using rs_mat = SparseMatrix<rational_type>;

        using lis_mat = SparseMatrix<lint_type>;
        using lrs_mat = SparseMatrix<lrational_type>;

        using is_vec = SparseVector<int_type>;
        using rs_vec = SparseVector<rational_type>;

        using lis_vec = SparseVector<lint_type>;
        using lrs_vec = SparseVector<lrational_type>;

    protected:
        using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;
        using LRatToLInt_Unary = Utility::RatToInt_Unary<lint_type, lrational_type>;

        RatToInt_Unary _rat_to_int;
        LRatToLInt_Unary _lrat_to_lint;

        static constexpr auto _to_single_int = conv<lint_type, int_type>::make;
        static constexpr auto _to_single_rat = conv<lrational_type, rational_type>::make;

    private:
        using my_type = base_alg<intT, rationalT, lintT, lrationalT>;

    public:
        logger &log;

        double tolerance = 0.000001;
        double err_tolerance = 0.01;
        size_t threads = 8;
        std::mt19937 *rgen = nullptr;

        base_alg(logger &a_log)
            : log(a_log), _rat_to_int(err_tolerance), _lrat_to_lint(err_tolerance)
        {
        }
    };
} // namespace ILPCuts
