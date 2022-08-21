#pragma once

#include "base_module.h"
#include "ntl_module.h"
#include <NTL/lzz_pX.h>
#include <NTL/vec_lzz_p.h>

namespace ILPCuts
{
    using namespace NTL;

    namespace ModuleNTL
    {
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

            vec_zz_p operator()(const vec_zz_p &b) const
            {
                size_t n = b.length();
                assert(n == _A.cols());

                is_vec sp_b = conv<vec_zz_p, is_vec>::make(b);

                sp_b = _A * sp_b;

                vec_zz_p res = conv<is_vec, vec_zz_p>::make(sp_b);

                return res;
            }
        };

        class WiedemannsAlg
        {
        private:
            vec_zz_p generate_nonzero_vec(size_t n)
            {
                size_t ntries = 0;
                int64_t p = zz_p::modulus();

                while (ntries < nonzero_lim)
                {
                    vec_zz_p res;
                    random(res, n);

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
            vec_zz_p generate_seq(size_t len, const oracleT &O, vec_zz_p b)
            {
                size_t n = b.length();

                vec_zz_p u = generate_nonzero_vec(n);
                // cout << "u = " << Utility::dense_vec_to_str(u) << endl;

                vec_zz_p res;
                res.SetLength(len);
                for (size_t i = 0; i < len; ++i)
                {
                    NTL::InnerProduct(res[i], u, b);
                    // res[i] = inner_product(u.cbegin(), u.cend(), b.cbegin(), zz_p(0));

                    b = O(b);
                }

                return res;
            }

            bool is_zero_vec(const vec_zz_p &x) const
            {
                return all_of(x.begin(), x.end(), [](const zz_p &e)
                              { return e == 0; });
            }

        public:
            logger &log;
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

            WiedemannsAlg(logger &alog)
                : log(alog) {}

            template <typename oracleT>
            vec_zz_p apply_poly(const oracleT &O, const zz_pX &fX, vec_zz_p b)
            {
                size_t n = b.length();
                size_t d = deg(fX);

                vec_zz_p res;
                res.SetLength(n, zz_p(0));
                for (size_t i = 0; i <= d; ++i)
                {
                    res += b * fX[i];
                    // transform(res.cbegin(), res.cend(), b.cbegin(), res.begin(), [&](const auto &e1, const auto &e2)
                    //           { return e1 + fX[i] * e2; });
                    b = O(b);
                }

                return res;
            }

            template <typename oracleT>
            zz_pX vec_min_poly(const oracleT &O, vec_zz_p b)
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
                    vec_zz_p seqX = generate_seq(len, O, b);

                    // out() << "seqX = " << seqX << endl;
                    // zz_pX seq_mX = seq_min_poly(seqX, len);
                    zz_pX seq_mX;
                    NTL::MinPolySeq(seq_mX, seqX, len / 2);

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
            vec_zz_p solve_kernel_v1(const oracleT &O)
            {
                size_t n = O.dim();
                int64_t p = zz_p::modulus();

                bool zero_root_poly_generated = false;
                zz_pX mX;
                vec_zz_p b;
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

    } // namespace ModuleNTL
} // namespace ILPCuts
