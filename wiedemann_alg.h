#pragma once

#include "base_module.h"
#include "ntl_module.h"
#include <NTL/lzz_pX.h>

namespace ILPCuts
{
    using namespace NTL;

    tuple<zz_pX, zz_pX> pade_approx(zz_pX bX, size_t n, size_t k)
    {
        zz_pX aX(INIT_MONO, n);

        zz_pX saX(INIT_MONO, 0, 1);
        zz_pX sbX;
        zz_pX taX;
        zz_pX tbX(INIT_MONO, 0, 1);

        while (deg(bX) != -1 && deg(bX) >= k)
        {
            zz_pX qX, rX;
            DivRem(qX, rX, aX, bX);

            zz_pX nsbX = saX - sbX * qX;
            zz_pX ntbX = taX - tbX * qX;

            saX = sbX;
            sbX = nsbX;

            taX = tbX;
            tbX = ntbX;

            aX = bX;
            bX = rX;

            bX.normalize();
        }

        tbX.normalize();

        return make_tuple(bX, tbX);
    }

    zz_pX seq_min_poly(const zz_pX &aX, size_t len)
    {
        size_t n = len / 2;

        zz_pX sX, tX;
        tie(sX, tX) = pade_approx(aX, 2 * n, n);

        size_t res_deg = max(1 + deg(sX), deg(tX));

        tX = reverse(tX);
        MakeMonic(tX);
        tX = LeftShift(tX, res_deg - deg(tX));

        return tX;
    }

    zz_pX seq_min_poly(const zz_pX &aX)
    {
        return seq_min_poly(aX, deg(aX) + 1);
    }

    template <typename intT>
    class Oracle_zzp
    {
    public:
        using int_type = intT;
        using is_mat = SparseMatrix<int_type>;

    private:
        using is_vec = SparseVector<int_type>;

        const is_mat _A;

    public:
        Oracle_zzp(is_mat A)
            : _A(move(A)) { assert(_A.cols() == _A.rows()); }

        size_t dim() const
        {
            return _A.cols();
        }

        vector<zz_p> operator()(const vector<zz_p> &b) const
        {
            size_t n = b.size();
            assert(n == _A.cols());

            is_vec sp_b = Utility::vec_conv<zz_p, int_type>(b);

            sp_b = _A * sp_b;

            vector<zz_p> res = Utility::vec_conv<int_type, zz_p>(sp_b);

            return res;
        }
    };

    class wiedemanns_alg
    {
    private:
        vector<zz_p> generate_nonzero_vec(size_t n)
        {
            size_t ntries = 0;
            int64_t p = zz_p::modulus();

            std::uniform_int_distribution<int64_t> ud(0, p - 1);

            while (ntries < nonzero_lim)
            {
                vector<zz_p> res(n, zz_p(0));

                for (size_t i = 0; i < n; ++i)
                {
                    res[i] = ud(*rgen);
                }

                if (!is_zero_vec(res))
                    return res;

                ++ntries;
            }

            stringstream msg;
            msg << "ILPCuts error: failed to generate nonzero vector in Wiedemann's algorithm:\n"
                << "limit: " << nonzero_lim << "\n"
                << "tries: " << ntries << "\n"
                << "modulus: " << p << "\n";

            throw Exception(msg.str());
        }

        template <typename oracleT>
        zz_pX generate_seq(size_t len, const oracleT &O, vector<zz_p> b)
        {
            size_t n = b.size();

            vector<zz_p> u = generate_nonzero_vec(n);
            // cout << "u = " << Utility::dense_vec_to_str(u) << endl;

            zz_pX res;
            res.SetLength(len);
            for (size_t i = 0; i < len; ++i)
            {

                // cout << "i = " << i << "; b = " << Utility::dense_vec_to_str(b) << endl;
                res[i] = inner_product(u.cbegin(), u.cend(), b.cbegin(), zz_p(0));
                // cout << "res[i] = " << res[i] << endl;

                b = O(b);
            }

            res.normalize();
            return res;
        }

        bool is_zero_vec(const vector<zz_p> &x) const
        {
            return all_of(x.begin(), x.end(), [](const zz_p &e)
                          { return e == 0; });
        }

    public:
        logger &log;
        std::mt19937 *rgen;
        size_t nonzero_lim = 3;
        size_t kernel_lim = 30;
        size_t vec_min_poly_lim = 30;
        bool silent = false;

        logger::stream &out()
        {
            if (!silent)
            {
                log.out() << "ILPCuts Wiedemann's alg.: ";
                return log.out();
            }

            return log.empty();
        }

        wiedemanns_alg(logger &alog)
            : log(alog) {}

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
            long p = zz_p::modulus();
            assert(n == b.size());

            zz_pX mX(INIT_MONO, 0, 1);

            size_t ntries = 0;
            while (!is_zero_vec(b) && ntries < vec_min_poly_lim)
            {
                // out() << "b = " << Utility::dense_vec_to_str(b) << endl;
                // out() << "mX = " << mX << "; deg = " << deg(mX) << endl;
                size_t len = 2 * (n - deg(mX));
                zz_pX seqX = generate_seq(len, O, b);

                // out() << "seqX = " << seqX << endl;
                zz_pX seq_mX = seq_min_poly(seqX, len);

                // out() << "seq_mX = " << seq_mX << endl;

                b = apply_poly(O, seq_mX, b);

                mX *= seq_mX;

                ++ntries;
            }

            if (!is_zero_vec(b))
            {
                stringstream msg;
                msg << "ILPCuts error: failed to compute minimal polynomial in Wiedemann's algorithm:\n"
                    << "limit: " << vec_min_poly_lim << "\n"
                    << "tries: " << ntries << "\n"
                    << "modulus: " << p << "\n";

                throw Exception(msg.str());
            }

            return mX;
        }

        template <typename oracleT>
        vector<zz_p> solve_kernel_v1(const oracleT &O)
        {
            size_t n = O.dim();
            int64_t p = zz_p::modulus();

            bool zero_root_poly_generated = false;
            zz_pX mX;
            vector<zz_p> b;
            size_t ntries = 0;
            while (!zero_root_poly_generated && ntries < kernel_lim)
            {
                b = generate_nonzero_vec(n);

                mX = vec_min_poly(O, b);

                if (mX[0] == 0)
                    zero_root_poly_generated = true;

                ++ntries;
            }

            if (!zero_root_poly_generated)
            {
                stringstream msg;
                msg << "ILPCuts error: failed to compute kernel vector in Wiedemann's algorithm:\n"
                    << "limit: " << kernel_lim << "\n"
                    << "tries: " << ntries << "\n"
                    << "modulus: " << p << "\n";

                throw Exception(msg.str());
            }

            mX = RightShift(mX, 1);
            mX.normalize();

            return apply_poly(O, mX, b);
        }
    };
} // namespace ILPCuts
