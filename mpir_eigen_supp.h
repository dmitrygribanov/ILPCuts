#pragma once

#include <mpirxx.h>

#include <Eigen/Core>

namespace Eigen {
    template<> struct NumTraits<mpq_class> : GenericNumTraits<mpq_class>
    {
        typedef mpq_class Real;
        typedef mpq_class NonInteger;
        typedef mpq_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 150,
            MulCost = 100
        };
    };

    template<> struct NumTraits<mpf_class> : GenericNumTraits<mpf_class>
    {
        typedef mpf_class Real;
        typedef mpf_class NonInteger;
        typedef mpf_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 75,
            MulCost = 50
        };
    };

    template<> struct NumTraits<mpz_class> : GenericNumTraits<mpz_class>
    {
        typedef mpz_class Real;
        typedef mpz_class NonInteger;
        typedef mpz_class Nested;

        static inline Real epsilon() { return 0; }
        static inline Real dummy_precision() { return 0; }
        static inline int digits10() { return 0; }

        enum {
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

mpq_class floor(const mpq_class& x) {
    mpz_class num = x.get_num();
    mpz_class den = x.get_den();
    
    return mpz_class(num/den);
}

mpq_class ceil(const mpq_class& x) {
    mpz_class num = x.get_num();
    mpz_class den = x.get_den();

    if ( num % den != 0 && x > 0 ) {
        return mpz_class(num / den + 1);
    }
    else {
        return mpz_class(num / den);
    }
}

mpf_class round(const mpf_class& x) {
    if (x > mpf_class(0))
        return floor(x + 0.5);
    else
        return ceil(x - 0.5);
    // return floor(x); // (x > mpf_class(0)) ? mpf_class(floor(x + mpf_class(0.5))) : mpf_class(ceil(x - mpf_class(0.5)));
}

mpq_class round(const mpq_class& x) {
    if (x >= mpq_class(0))
        return floor(mpq_class(x + mpq_class(1, 2)));
    else
    {
        //cout << "x: " << x << endl;
        return ceil(mpq_class(x - mpq_class(1, 2)));
    }
    //return mpq_class(ceil(x));
    // return floor(x); // (x > mpf_class(0)) ? mpf_class(floor(x + mpf_class(0.5))) : mpf_class(ceil(x - mpf_class(0.5)));
}

