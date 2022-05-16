#pragma once

#include "conic_cuts_base.h"
#include "linbox_module.h"

namespace ILPCuts
{
    template <typename intT, typename rationalT, typename lintT, typename lratT>
    class ConicCutsGenerator_Zp : public ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>
    {
    private:
        using base = ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>;

    public:
        using typename base::int_type;
        using typename base::lint_type;
        using typename base::lrational_type;
        using typename base::rational_type;

        using typename base::is_mat;
        using typename base::lis_mat;
        using typename base::lrs_mat;
        using typename base::rs_mat;

        using typename base::is_vec;
        using typename base::lis_vec;
        using typename base::lrs_vec;
        using typename base::rs_vec;

        using typename base::cone_traits;
        using typename base::dir_coeffs;
        using typename base::lsparce_transform;
        using typename base::process_func;
        using typename base::sparce_transform;
        using typename base::transform_base;

        using typename base::base;
        using typename base::cone_traits;

        size_t max_conservative_prime = 11;

    private:
        std::vector<int_type> _primes_data;

        std::unique_ptr<cone_traits> create_cone_traits(const is_mat &A, const lint_type &det, const std::vector<rs_vec> &directions_pool)
        {

            auto res = std::make_unique<cone_traits>();

            res->A = A;
            res->det = det;

            SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU(A.template cast<rational_type>());

            for (size_t i = 0; i < directions_pool.size(); ++i)
            {
                rs_vec x = decomp_LU.solve(directions_pool[i]);
                res->dir_coeffs_pull.emplace_back(dir_coeffs{x, i});
            }

            return std::move(res);
        }

        void process_odd_conservative_case_zp(const cone_traits &C, process_func pf, int_type p)
        {
            is_vec x = Utility::LinBox::compute_some_kernel_vec(C.A, p);

            for (typename is_vec::InnerIterator it(x); it; ++it)
            {
                if (it.value() > 0)
                    pf(sparce_transform(this, static_cast<size_t>(it.index()), std::move(x), p));
            }
        }

        void process_even_case(const cone_traits &C, process_func pf)
        {
            process_odd_conservative_case_zp(C, pf, 2);
        }

        void process_odd_case_zp(const cone_traits &C, process_func pf, int_type p)
        {
            is_vec x = Utility::LinBox::compute_some_kernel_vec(C.A, p);

            for (typename lis_vec::InnerIterator it(x); it; ++it)
            {
                if (it.value() > 0)
                {
                    int_type lambda = Utility::mod_equation_solve<int_type>(it.value(), p - gcd(p, it.value()), p);
                    is_vec t = (x * lambda).unaryExpr(Utility::mod_unaryop(p));

                    pf(sparce_transform(this, static_cast<size_t>(it.index()), std::move(t), p));
                }
            }
        }

        void process_odd_case_basic(const cone_traits &C, process_func pf) {}

        void process_node(const cone_traits &C, process_func pf) override
        {
            int_type p = Utility::primes::find_divisor(C.det, _primes_data);

            if (p == 2)
            {
                out() << "even step..." << std::endl;
                process_even_case(C, pf);
            }
            else if (p <= max_conservative_prime)
            {
                out() << "odd conservative zp step..." << std::endl;
                process_odd_conservative_case_zp(C, pf, p);
            }
            else if (p != 0)
            {
                out() << "odd zp step..." << std::endl;
                process_odd_case_zp(C, pf, p);
            }
            else
            {
                out() << "odd basic step..." << std::endl;
                process_odd_case_basic(C, pf);
            }
        }

    public:
        std::vector<is_vec> generate_cuts(const is_mat &A, const rs_vec &optLP, const lint_type &det)
        {
            auto directions_pool = std::move(base::generate_directions_pool(A));

            size_t n = A.cols();
            size_t primes_num = std::min(Utility::primes::max_download, n);
            _primes_data = Utility::primes::download<int_type>(primes_num);

            return std::move(base::generate_cuts(
                create_cone_traits(A, det, directions_pool),
                optLP));
        }

        std::vector<is_vec> generate_cuts(const is_mat &A, const rs_vec &optLP) override
        {
            try
            {
                auto start = std::chrono::steady_clock::now();
                lint_type det = Utility::LinBox::compute_det(A.template cast<lint_type>());
                auto end = std::chrono::steady_clock::now();

                std::chrono::duration<double> elapsed_seconds = end - start;

                out() << "ConicCuts: determinant computation time: " << elapsed_seconds.count() << "s\n";

                return generate_cuts(A, optLP, det);
            }
            catch (ToleranceErrorException exc)
            {
                out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
                err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
            }
        }
    };
}