#pragma once

#include <eigen3/Eigen/SparseCore>

#include "gmp_supp.h"

namespace ILPCuts
{
    using namespace Eigen;
    using namespace std;

    using int_type = int64_t;
    using lint_type = mpz_class;
    using rational_type = double;
    using lrational_type = mpq_class;

    class Exception : public std::runtime_error
    {
    public:
        Exception()
            : std::runtime_error("ILPCuts lib error") {}

        Exception(std::string msg)
            : std::runtime_error(msg) {}
    };

    class ToleranceErrorException : public Exception
    {
    public:
        ToleranceErrorException()
            : Exception("ILPCuts lib error: a tolerance-error has occurred.") {}

        ToleranceErrorException(std::string msg)
            : Exception(msg) {}
    };

    template <typename inT, typename outT>
    struct conv
    {
        static outT make(const inT &a)
        {
            return outT(a);
        }
    };

    template <typename inT, typename outVecT>
    struct conv<SparseVector<inT>, outVecT>
    {
        static outVecT make(const SparseVector<inT> &a)
        {
            using out_val_type = typename outVecT::value_type;

            size_t n = a.size();

            vector<out_val_type> res(n, out_val_type(0));

            for (typename SparseVector<inT>::InnerIterator it(a); it; ++it)
            {
                res[it.index()] = conv<inT, out_val_type>::make(it.value());
            }

            return res;
        }
    };

    template <typename inVecT, typename outT>
    struct conv<inVecT, SparseVector<outT>>
    {
        static SparseVector<outT> make(const inVecT &a)
        {
            using in_val_type = typename inVecT::value_type;

            size_t n = a.size();

            SparseVector<outT> res(n);
            for (size_t i = 0; i < n; ++i)
            {
                if (a[i] != 0)
                    res.insertBack(i) = conv<in_val_type, outT>::make(a[i]);
            }

            return res;
        }
    };

    namespace Utility
    {
        template <typename intT, typename ratT>
        intT rat_to_int(const ratT &rat_x, double err_tolerance)
        {
            ratT rounded_x = round(rat_x);
            ratT diff = abs(rat_x - rounded_x);
            if (diff <= err_tolerance)
            {
                return conv<ratT, intT>::make(rounded_x);
            }
            else
            {
                std::stringstream msg;
                msg << "ILPCuts tolerance error:\n diff = (" << diff
                    << ")\n value = (" << rat_x
                    << ")\n rounded value = (" << rounded_x << ")";
                throw ToleranceErrorException(msg.str());
            }
        }

        template <typename intT>
        intT mod2_unaryop(const intT &x)
        {
            return abs(x % intT(2));
        }

        template <>
        int64_t mod2_unaryop(const int64_t &x)
        {
            return x & 1;
        }

        template <typename intT>
        class mod_unaryop
        {
        private:
            const intT &_mod;

        public:
            mod_unaryop(const intT &mod)
                : _mod(mod) {}

            intT operator()(const intT &x) const
            {
                intT res(x % _mod);
                if (res < 0)
                    res += _mod;
                return res;
            }
        };

        template <typename ratT>
        class Round_Unary
        {
        public:
            ratT operator()(const ratT &rat_x) const
            {
                return round(rat_x);
            }
        };

        template <typename ratT>
        class Floor_Unary
        {
        public:
            ratT operator()(const ratT &rat_x) const
            {
                return floor(rat_x);
            }
        };

        template <typename ratT>
        class Ceil_Unary
        {
        public:
            ratT operator()(const ratT &rat_x) const
            {
                return ceil(rat_x);
            }
        };

        template <typename intT, typename ratT>
        class RatToInt_Unary
        {
        private:
            const double &_err_tolerance;

        public:
            explicit RatToInt_Unary(const double &err_tolerance)
                : _err_tolerance(err_tolerance) {}

            intT operator()(const ratT &rat_x) const
            {
                return rat_to_int<intT, ratT>(rat_x, _err_tolerance);
            }

            const double &get_err_tolerance() const { return _err_tolerance; }
        };

        template <typename T>
        T euclid_bezout(
            const T &a,
            const T &b,
            T &u,
            T &v)
        {
            u = 1;
            v = 0;
            T u_new = v;
            T v_new = u;
            T r = abs(a);
            T r_new = abs(b);

            while (r_new != 0)
            {
                T q = r / r_new;

                T aux = u_new;
                u_new = u - q * u_new;
                u = aux;

                aux = v_new;
                v_new = v - q * v_new;
                v = aux;

                aux = r_new;
                r_new = r - q * r_new;
                r = aux;
            }

            return r;
        }

        template <typename intT>
        intT mod_equation_solve(const intT &a, const intT &b, const intT &mod);

        template <>
        int64_t mod_equation_solve(const int64_t &a, const int64_t &b, const int64_t &mod)
        {
#ifdef _DEBUG
            assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif

            int64_t x = 0, y = 0;
            int64_t gcd_a_mod = euclid_bezout(a, mod, x, y);

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

    using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;
    using LRatToLInt_Unary = Utility::RatToInt_Unary<lint_type, lrational_type>;

    constexpr auto to_single_int = conv<lint_type, int_type>::make;
    constexpr auto to_single_rat = conv<lrational_type, rational_type>::make;

    using is_mat = SparseMatrix<int_type>;
    using rs_mat = SparseMatrix<rational_type>;

    using lis_mat = SparseMatrix<lint_type>;
    using lrs_mat = SparseMatrix<lrational_type>;

    using is_vec = SparseVector<int_type>;
    using rs_vec = SparseVector<rational_type>;

    using lis_vec = SparseVector<lint_type>;
    using lrs_vec = SparseVector<lrational_type>;

} // namespace ILPCuts
