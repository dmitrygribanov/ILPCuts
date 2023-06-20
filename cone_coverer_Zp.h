#pragma once

#include "cone_coverer_base.h"
#include "ntl_module.h"
#include "primes_downloader.h"

namespace ILPCuts
{
    class ConeCoverer_Zp : public ConeCoverer_Base
    {
    private:
        using base = ConeCoverer_Base;

    public:
        size_t max_conservative_prime = 101;

    private:
        std::vector<int_type> _primes_data;

        std::unique_ptr<cone_traits> create_cone_traits(const is_mat &A, const lint_type &det, const std::vector<rs_vec> &directions_pool)
        {

            auto res = std::make_unique<cone_traits>();

            res->A = A;
            res->det = det;

            SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU(A.cast<rational_type>());

            for (size_t i = 0; i < directions_pool.size(); ++i)
            {
                rs_vec x = decomp_LU.solve(directions_pool[i]);
                res->dir_coeffs_pull.emplace_back(dir_coeffs{x, i});
            }

            return std::move(res);
        }

        void process_odd_conservative_case_zp(const cone_traits &C, process_func pf, int_type p)
        {
            is_vec x = ModuleNTL::kernel_vec_Zp(C.A, p);

            for (typename is_vec::InnerIterator it(x); it; ++it)
            {
                if (it.value() > 0)
                    pf(sparse_transform(tolerance, err_tolerance, static_cast<size_t>(it.index()), std::move(x), p));
            }
        }

        void process_even_case(const cone_traits &C, process_func pf)
        {
            is_vec x = ModuleNTL::kernel_vec_GF2(C.A);

            for (typename is_vec::InnerIterator it(x); it; ++it)
            {
                if (it.value() > 0)
                    pf(sparse_transform(tolerance, err_tolerance, static_cast<size_t>(it.index()), std::move(x), 2));
            }
        }

        template <typename int_typeT, typename sparse_transformT>
        void process_odd_x(const SparseVector<int_typeT> &x, process_func pf, const int_typeT &denum)
        {
            using vec_type = SparseVector<int_typeT>;

            for (typename vec_type::InnerIterator it(x); it; ++it)
            {
                if (it.value() > 0)
                {
                    int_typeT lambda = Utility::mod_equation_solve<int_typeT>(it.value(), denum - gcd(denum, it.value()), denum);
                    vec_type t = (x * lambda).unaryExpr(Utility::mod_unaryop(denum));

                    pf(sparse_transformT(tolerance, err_tolerance, static_cast<size_t>(it.index()), std::move(t), denum));
                }
            }
        }

        void process_odd_case_zp(const cone_traits &C, process_func pf, int_type p)
        {
            is_vec x = ModuleNTL::kernel_vec_Zp(C.A, p);

            process_odd_x<int_type, sparse_transform>(x, pf, p);
        }

        void process_odd_case_basic(const cone_traits &C, process_func pf)
        {
            const lint_type &det = C.det;
            size_t n = C.A.cols();

            bool has_non_zero_solution = false;
            lis_vec x;
            while (!has_non_zero_solution)
            {
                lis_vec b = Utility::bernouly_vec_generate<int_type>(n, *base::rgen).cast<lint_type>();

                x = ModuleNTL::solve_ZZ(C.A.cast<lint_type>().eval(), b);

                x = x.unaryExpr(Utility::mod_unaryop(det)).eval();
                Utility::prune_zeroes_int(x);

                if (x.nonZeros() > 0)
                    has_non_zero_solution = true;
            }

            process_odd_x<lint_type, lsparse_transform>(x, pf, det);
        }

        void process_node(const cone_traits &C, process_func pf) override
        {
            int_type p = Utility::primes::find_divisor(C.det, _primes_data);

            if (p == 2)
            {
                out() << "even step...\n";
                process_even_case(C, pf);
            }
            else if (p != 0 && p <= max_conservative_prime)
            {
                out() << "odd conservative zp step...\n";
                process_odd_conservative_case_zp(C, pf, p);
            }
            else if (p != 0)
            {
                out() << "odd zp step...\n";
                process_odd_case_zp(C, pf, p);
            }
            else
            {
                out() << "odd basic step...\n";
                process_odd_case_basic(C, pf);
            }
        }

    public:
        using base::base;

        std::vector<is_mat> generate_cuts(const is_mat &A, const lint_type &det)
        {
            out() << "Generating directions pool ..." << endl;

            auto directions_pool = std::move(base::generate_directions_pool(A));

            size_t n = A.cols();
            Timer timer;
            timer.start();

            out() << "Primes downloading is started ..." << endl;

            size_t primes_num = std::min(Utility::primes::max_download, n * n);
            _primes_data = Utility::primes::download<int_type>(primes_num);

            out() << "Primes downloading is finished, time = " << timer.finish_str();

            out() << "Constructing initial cone ..." << endl;
            timer.start();

            auto initial_cone = create_cone_traits(A, det, directions_pool);

            out() << "Initial cone construction is finished, time = " << timer.finish_str();

            return base::generate_cuts(std::move(initial_cone));
        }

        std::vector<is_mat> generate_cuts(const is_mat &A) override
        {
            try
            {
                Timer timer;
                timer.start();

                out() << "Determinant computation is started ..." << endl;
                // lint_type det = _lrat_to_lint(Utility::compute_abs_det(A.cast<lrational_type>().eval()));
                lint_type det = ModuleNTL::det(A.cast<lint_type>().eval());
                // lint_type det;
                // ModuleFlint::compute_det(A.cast<lint_type>().eval(), det);
                // //  lint_type det = Utility::LinBox::compute_det(A.cast<lint_type>().eval());
                out() << "determinant computation time: " << timer.finish_str()
                      << "determinant is: " << det << endl;

                return generate_cuts(A, abs(det));
            }
            catch (ToleranceErrorException exc)
            {
                out() << "determinant computing tolerance error: " << exc.what() << endl;
                return std::vector<is_mat>();
            }
        }
    };
} // namespace ILPCuts