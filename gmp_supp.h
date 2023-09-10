#pragma once

#include <gmpxx.h>

#include "dependencies/Eigen/Core"
namespace Eigen
{
    template <>
    struct NumTraits<mpq_class> : GenericNumTraits<mpq_class>
    {
        typedef mpq_class Real;
        typedef mpq_class NonInteger;
        typedef mpq_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 150,
            MulCost = 100
        };
    };

    template <>
    struct NumTraits<mpf_class> : GenericNumTraits<mpf_class>
    {
        typedef mpf_class Real;
        typedef mpf_class NonInteger;
        typedef mpf_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum
        {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 50,
            MulCost = 75
        };
    };

    template <>
    struct NumTraits<mpz_class> : GenericNumTraits<mpz_class>
    {
        typedef mpz_class Real;
        typedef mpz_class NonInteger;
        typedef mpz_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

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

mpq_class floor(const mpq_class &x)
{
    mpz_class num = x.get_num();
    mpz_class den = x.get_den();

    return mpz_class(num / den);
}

mpq_class ceil(const mpq_class &x)
{
    mpz_class num = x.get_num();
    mpz_class den = x.get_den();

    if (num % den != 0 && x > 0)
    {
        return mpz_class(num / den + 1);
    }
    else
    {
        return mpz_class(num / den);
    }
}

mpf_class round(const mpf_class &x)
{
    if (x > mpf_class(0))
        return floor(x + 0.5);
    else
        return ceil(x - 0.5);
    // return floor(x); // (x > mpf_class(0)) ? mpf_class(floor(x + mpf_class(0.5))) : mpf_class(ceil(x - mpf_class(0.5)));
}

mpq_class round(const mpq_class &x)
{
    if (x >= mpq_class(0))
        return floor(x + mpq_class(1, 2));
    else
    {
        return ceil(x - mpq_class(1, 2));
    }
    // return mpq_class(ceil(x));
    //  return floor(x); // (x > mpf_class(0)) ? mpf_class(floor(x + mpf_class(0.5))) : mpf_class(ceil(x - mpf_class(0.5)));
}

namespace ILPCuts
{
    template <typename inT, typename outT>
    struct conv;
    template <>
    struct conv<mpz_class, int64_t>
    {
        static int64_t make(const mpz_class &a)
        {
            return a.get_si();
        }
    };

    template <>
    struct conv<mpf_class, double>
    {
        static double make(const mpf_class &a)
        {
            return a.get_si();
        }
    };

    template <>
    struct conv<mpq_class, double>
    {
        static double make(const mpq_class &a)
        {
            return a.get_d();
        }
    };

    template <>
    struct conv<mpf_class, int64_t>
    {
        static int64_t make(const mpf_class &a)
        {
            return static_cast<int64_t>(a.get_d());
        }
    };

    template <>
    struct conv<mpq_class, int64_t>
    {
        static int64_t make(const mpq_class &a)
        {
            return mpz_class(a.get_num() / a.get_den()).get_si();
        }
    };

    namespace Utility
    {
        template <typename intT>
        intT mod_equation_solve(const intT &a, const intT &b, const intT &mod);

        template <>
        mpz_class mod_equation_solve(const mpz_class &a, const mpz_class &b, const mpz_class &mod)
        {
#ifdef _DEBUG
            assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif

            mpz_class x = 0, y = 0;
            mpz_class gcd_a_mod = 0; // = Internal::euclid_bezout(a, mod, x, y);

            mpz_gcdext(gcd_a_mod.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

            if (x < 0)
                x += mod;

#ifdef _DEBUG
            assert(x > 0);
            assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
            x *= (b / gcd_a_mod);

            return x;
        }

    } // namespace Utility

} // namespace ILPCuts