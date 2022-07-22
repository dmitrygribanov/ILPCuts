#pragma once

#include "boost_tee_logging.h"

namespace ILPCuts
{
    using namespace std;

    template <typename intT, typename lintT, typename ratT, typename lratT>
    void print_numeric_types_info(logger &log)
    {
        log.out() << "ILPCuts numeric-types info:\n"
                  << "short-int type: " << typeid(intT).name() << "\n"
                  << "long-int type: " << typeid(lintT).name() << "\n"
                  << "short-rational type: " << typeid(ratT).name() << "\n"
                  << "long-rational type: " << typeid(lratT).name() << "\n"
                  << endl;
    }

}