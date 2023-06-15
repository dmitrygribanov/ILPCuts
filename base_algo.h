#pragma once

#include <numeric>
#include <functional>
#include <memory>

#include "basic_utils.h"
#include "gmp_supp.h"
#include "timer.h"
#include "sparse_transform.h"
#include "boost_tee_logging.h"

namespace ILPCuts
{
#define SINGLE_NUMERIC_TYPES_DEFS(INT_TYPE, RATIONAL_TYPE) \
    using rational_type = RATIONAL_TYPE;                   \
    using int_type = INT_TYPE;                             \
                                                           \
    using is_mat = SparseMatrix<int_type>;                 \
    using rs_mat = SparseMatrix<rational_type>;            \
                                                           \
    using is_vec = SparseVector<int_type>;                 \
    using rs_vec = SparseVector<rational_type>;

#define LONG_NUMERIC_TYPES_DEFS(LINT_TYPE, LRATIONAL_TYPE) \
    using lrational_type = LRATIONAL_TYPE;                 \
    using lint_type = LINT_TYPE;                           \
                                                           \
    using lis_mat = SparseMatrix<lint_type>;               \
    using lrs_mat = SparseMatrix<lrational_type>;          \
                                                           \
    using lis_vec = SparseVector<lint_type>;               \
    using lrs_vec = SparseVector<lrational_type>;
    template <typename intT, typename rationalT>
    class single_type_base_alg
    {
    public:
        SINGLE_NUMERIC_TYPES_DEFS(intT, rationalT)

    protected:
        using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;

        RatToInt_Unary _rat_to_int;

    private:
        using my_type = single_type_base_alg<intT, rationalT>;

    public:
        logger &log;

        double tolerance = 0.000001;
        double err_tolerance = 0.01;
        size_t threads = 8;
        std::mt19937 *rgen = nullptr;

        single_type_base_alg(logger &a_log)
            : log(a_log), _rat_to_int(err_tolerance)
        {
        }
    };
    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class dual_type_base_alg
    {
    public:
        SINGLE_NUMERIC_TYPES_DEFS(intT, rationalT)
        LONG_NUMERIC_TYPES_DEFS(lintT, lrationalT)

    protected:
        using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;
        using LRatToLInt_Unary = Utility::RatToInt_Unary<lint_type, lrational_type>;

        RatToInt_Unary _rat_to_int;
        LRatToLInt_Unary _lrat_to_lint;

        static constexpr auto _to_single_int = conv<lint_type, int_type>::make;
        static constexpr auto _to_single_rat = conv<lrational_type, rational_type>::make;

    private:
        using my_type = dual_type_base_alg<intT, rationalT, lintT, lrationalT>;

    public:
        logger &log;

        double tolerance = 0.000001;
        double err_tolerance = 0.01;
        size_t threads = 8;
        std::mt19937 *rgen = nullptr;

        dual_type_base_alg(logger &a_log)
            : log(a_log), _rat_to_int(err_tolerance), _lrat_to_lint(err_tolerance)
        {
        }
    };

#define SINGLE_TYPE_ALGO_DEFS           \
    using typename base::int_type;      \
    using typename base::rational_type; \
                                        \
    using typename base::is_mat;        \
    using typename base::rs_mat;        \
                                        \
    using typename base::is_vec;        \
    using typename base::rs_vec;        \
                                        \
    using base::_rat_to_int;            \
    using base::base;                   \
    using base::err_tolerance;          \
    using base::log;                    \
    using base::rgen;                   \
    using base::threads;                \
    using base::tolerance;

#define DUAL_TYPE_ALGO_DEFS              \
    SINGLE_TYPE_ALGO_DEFS                \
    using typename base::lint_type;      \
    using typename base::lrational_type; \
                                         \
    using typename base::lis_mat;        \
    using typename base::lrs_mat;        \
                                         \
    using typename base::lis_vec;        \
    using typename base::lrs_vec;        \
                                         \
    using base::_lrat_to_lint;           \
    using base::_to_single_int;          \
    using base::_to_single_rat;
} // namespace ILPCuts
