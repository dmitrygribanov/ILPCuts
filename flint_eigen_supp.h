#pragma once

#include <flint/fmpzxx.h>
#include <flint/fmpqxx.h>

#include <eigen3/Eigen/Core>

namespace Eigen
{
    template <>
    struct NumTraits<flint::fmpqxx> : GenericNumTraits<flint::fmpqxx>
    {
        typedef flint::fmpqxx Real;
        typedef flint::fmpqxx NonInteger;
        typedef flint::fmpqxx Nested;

        static inline Real epsilon() { return flint::fmpqxx(0l, 1ul); }
        static inline Real dummy_precision() { return flint::fmpqxx(0l, 1ul); }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 10,
            AddCost = 150,
            MulCost = 75
        };
    };

    // template <>
    // struct NumTraits<mpf_class> : GenericNumTraits<mpf_class>
    // {
    //     typedef mpf_class Real;
    //     typedef mpf_class NonInteger;
    //     typedef mpf_class Nested;

    //     static inline Real epsilon() { return 0; }
    //     static inline Real dummy_precision() { return 0; }
    //     static inline int digits10() { return 0; }

    //     enum
    //     {
    //         IsInteger = 0,
    //         IsSigned = 1,
    //         IsComplex = 0,
    //         RequireInitialization = 1,
    //         ReadCost = 6,
    //         AddCost = 75,
    //         MulCost = 50
    //     };
    // };

    template <>
    struct NumTraits<flint::fmpzxx> : GenericNumTraits<flint::fmpzxx>
    {
        typedef flint::fmpzxx Real;
        typedef flint::fmpzxx NonInteger;
        typedef flint::fmpzxx Nested;

        static inline Real epsilon() { return Real(0); }
        static inline Real dummy_precision() { return Real(0); }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 1,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 75,
            MulCost = 50
        };
    };
}