#pragma once

#include <vector>

// #include "numeric_module.h"
#include "Eigen_solving_utilities.h"
#include "ntl_module.h"

namespace ILPCuts
{
    class local_instance
    {
    public:
        enum class var_ILP_type
        {
            Int,
            Rat
        };

        using decomp_type = SparseLU<rs_mat, COLAMDOrdering<int>>;

        std::unique_ptr<lint_type> det;
        std::unique_ptr<decomp_type> decomp_LU;
        // logger *log;

        local_instance(/*logger *alog = nullptr*/) /*: log(alog)*/ {}

        local_instance(/*logger *alog,*/ size_t n)
            : /*log(alog),*/ A(n, n), c(n), v(n), varsLB(n), varsUB(n), ILP_var_types(n)
        {
        }

        local_instance(local_instance &&u)
            : /*log(alog),*/ A(std::move(u.A)), c(std::move(u.c)), v(std::move(u.v)),
              ILP_var_types(std::move(u.ILP_var_types)), varsLB(std::move(u.varsLB)), varsUB(std::move(u.varsUB)),
              det(std::move(u.det)), decomp_LU(std::move(u.decomp_LU))
        {
        }

        local_instance &operator=(local_instance &&u)
        {
            // log = u.log;
            A = std::move(u.A);
            c = std::move(u.c);
            v = std::move(u.v);
            ILP_var_types = std::move(u.ILP_var_types);
            varsLB = std::move(u.varsLB);
            varsUB = std::move(u.varsUB);
            det = std::move(u.det);
            decomp_LU = std::move(u.decomp_LU);

            return *this;
        }

        rs_mat A;
        // is_vec b;
        rs_vec c;
        rs_vec v;
        std::vector<var_ILP_type> ILP_var_types;
        rs_vec varsLB;
        rs_vec varsUB;

        size_t get_num_vars() const { return A.cols(); }

        // logger::stream &out()
        // {
        //     log->out() << "ILPCuts-LocalInstance: ";
        //     return log->out();
        // }

        const lint_type &get_det(logger &log, double err_tolerance)
        {
            if (!det)
            {
                LRatToLInt_Unary lrat_to_lint(err_tolerance);
                // check_decomposition(log);

                log.out() << "ILPCuts-LocalInstance: determinant computation is started ..." << endl;
                Timer timer;
                timer.start();

                // det = std::make_unique<lint_type>(lrat_to_lint(decomp_LU->absDeterminant()));
                //  lint_type det = _lrat_to_lint(Utility::compute_abs_det(A.cast<lrational_type>().eval()));
                det = std::make_unique<lint_type>(ModuleNTL::det(A.cast<lint_type>().eval()));
                // lint_type det;
                // ModuleFlint::compute_det(A.cast<lint_type>().eval(), det);
                // //  lint_type det = Utility::LinBox::compute_det(A.cast<lint_type>().eval());

                log.out() << "determinant computation time: " << timer.finish_str()
                          << "determinant is: " << *det << endl;
            }

            return *det;
        }

        void check_decomposition(logger &log)
        {
            if (!decomp_LU)
            {
                Timer timer;
                timer.start();
                log.out() << "ILPCuts-LocalInstance: LU decomposition computation is started..." << endl;

                decomp_LU = std::make_unique<decomp_type>(A);

                if (decomp_LU->info() != ComputationInfo::Success)
                    throw Utility::EigenSolveException(decomp_LU->info());

                log.out() << "LU decomposition computation time: " << timer.finish_str();
            }
        }

        rs_mat solve_AX_eq_B(logger &log, const rs_mat &B)
        {
            check_decomposition(log);

            return decomp_LU->solve(B);
        }

        rs_mat solve_ATX_eq_B(logger &log, const rs_mat &B)
        {
            check_decomposition(log);

            return decomp_LU->transpose().solve(B);
        }

        rs_vec solve_AX_eq_b(logger &log, const rs_vec &b)
        {
            check_decomposition(log);

            return decomp_LU->solve(b);
        }

        rs_vec solve_ATX_eq_b(logger &log, const rs_vec &b)
        {
            check_decomposition(log);

            return decomp_LU->transpose().solve(b);
        }

        void adjust_normal_cone(logger &log, const rs_vec &c, double tolerance)
        {
            check_decomposition(log);

            bool is_updated = false;

            std::tie(is_updated, A) = Utility::adjust_normal_cone(std::move(A), c, tolerance, *decomp_LU);

            if (is_updated)
            {
                decomp_LU.release();
            }
        }

        void adjust_normal_cone(logger &log, double tolerance)
        {
            adjust_normal_cone(log, c, tolerance);
        }

        // size_t get_num_constr() const { return A.rows(); }
    };
}