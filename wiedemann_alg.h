#pragma once

#include "base_module.h"
#include "ntl_module.h"
#include <NTL/lzz_pX.h>

namespace ILPCuts
{
    using namespace NTL;

    tuple<zz_pX, zz_pX> pade_approx(zz_pX bX, size_t k)
    {

        size_t n = deg(bX) + 1;
        zz_pX aX(INIT_MONO, n);

        zz_pX saX(INIT_MONO, 0, 1);
        zz_pX sbX;
        zz_pX taX;
        zz_pX tbX(INIT_MONO, 0, 1);

        while (deg(bX) >= k)
        {
            zz_pX qX, rX;
            DivRem(qX, rX, aX, bX);

            zz_pX nsaX = saX - sbX * qX;
            zz_pX ntbX = taX - tbX * qX;

            saX = sbX;
            sbX = nsaX;

            taX = tbX;
            tbX = ntbX;

            aX = bX;
            bX = rX;
        }

        return make_tuple(bX, tbX);
    }

    zz_pX seq_min_poly(const zz_pX &aX)
    {
        size_t n = (deg(aX) + 1) / 2;

        zz_pX sX, tX;
        tie(sX, tX) = pade_approx(aX, n);

        size_t res_deg = max(1 + deg(sX), deg(tX));

        tX = reverse(tX);
        MakeMonic(tX);
        tX = LeftShift(tX, res_deg - deg(tX));

        return tX;
    }

    template <typename intT>
    class Oracle_zzp
    {
    public:
        using int_type = intT;
        using is_mat = SparseMatrix<int_type>;

    private:
        using is_vec = SparseVector<int_type>;

        const is_mat &_A;

    public:
        Oracle_zzp(const is_mat &A)
            : _A(A) { assert(_A.cols() == _A.rows()); }

        size_t dim() const
        {
            return _A.cols();
        }

        vector<zz_p> operator()(const vector<zz_p> &b) const
        {
            size_t n = b.size();
            assert(n == _A.cols());

            is_vec sp_b(n);
            for (size_t i = 0; i < n; ++i)
            {
                if (b[i] != 0)
                    sp_b.insertBack(i) = rep(b[i]);
            }
            sp_b = _A * sp_b;

            vector<zz_p> res(n, zz_p(0));
            for (typename is_vec::InnerIterator it(sp_b); it; ++it)
            {
                res[it.index()] = it.value();
            }

            return res;
        }
    };

    class wiedemanns_alg
    {
    private:
        std::mt19937 &_rgen;

        template <typename oracleT>
        zz_pX generate_seq(size_t len, const oracleT &O, vector<zz_p> b)
        {
            size_t n = b.size();

            vector<zz_p> u;
            generate_n(back_inserter(u), n, random_zz_p);

            zz_pX res;
            res.SetLength(len);
            for (size_t i = 0; i < len; ++i)
            {
                res[i] = inner_product(u.begin(), u.end(), b.begin(), zz_p(0));

                b = O(b);
            }

            return res;
        }

        bool is_zero_vec(const vector<zz_p> &x) const
        {
            return all_of(x.begin(), x.end(), [](const zz_p &e)
                          { return e == 0; });
        }

    public:
        logger &log;

        wiedemanns_alg(std::mt19937 &rgen, logger &alog)
            : _rgen(rgen), log(alog) {}

        template <typename oracleT>
        vector<zz_p> apply_poly(const oracleT &O, const zz_pX &fX, vector<zz_p> b)
        {
            size_t n = b.size();
            size_t d = deg(fX);

            vector<zz_p> res(n, zz_p(0));
            for (size_t i = 0; i <= d; ++i)
            {
                transform(res.cbegin(), res.cend(), b.cbegin(), res.begin(), [&](const auto &e1, const auto &e2)
                          { return e1 + fX[i] * e2; });
                b = O(b);
            }

            return res;
        }

        template <typename oracleT>
        zz_pX vec_min_poly(const oracleT &O, vector<zz_p> b)
        {
            size_t n = O.dim();
            assert(n == b.size());

            long p = zz_p::modulus();

            zz_pX mX(INIT_MONO, 0, 1);

            while (is_zero_vec(b))
            {
                zz_pX seqX = generate_seq(2 * (n - deg(mX)), O, b);

                zz_pX seq_mX = seq_min_poly(seqX);

                b = apply_poly(O, seq_mX, b);

                mX *= seq_mX;
            }

            return mX;
        }

        template <typename oracleT>
        vector<zz_p> solve_kernel(const oracleT &O)
        {
            size_t n = O.dim();

            bool zero_root_poly_generated = false;
            zz_pX mX(INIT_MONO, 0, 1);
            vector<zz_p> b;
            while (!zero_root_poly_generated)
            {
                vector<zz_p> x;
                generate_n(back_inserter(x), n, random_zz_p);

                if (is_zero_vec(x))
                    continue;

                b = O(x);
                if (is_zero_vec(b))
                    return x;

                mX = vec_min_poly(O, b);
            }
        }
    };
} // namespace ILPCuts
