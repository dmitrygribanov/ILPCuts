#pragma once

#include <numeric>
#include <functional>
#include <memory>

#include "timer.h"
#include "basic_utils.h"
#include "numeric_module.h"
#include "boost_tee_logging.h"

namespace ILPCuts
{
    class base_alg
    {
    public:
        logger &log;

        RatToInt_Unary _rat_to_int;
        LRatToLInt_Unary _lrat_to_lint;

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
