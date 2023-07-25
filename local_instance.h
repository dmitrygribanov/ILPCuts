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
            LRatToLInt_Unary lrat_to_lint(err_tolerance);

            Timer timer;
            timer.start();

            if (!det)
            {
                if (!decomp_LU)
                {
                    decomp_LU = std::make_unique<decomp_type>(A);
                }

                det = std::make_unique<lint_type>(lrat_to_lint(decomp_LU->absDeterminant()));
            }

            return *det;
            // out() << "Determinant computation is started ..." << endl;

            // lint_type det = _lrat_to_lint(Utility::compute_abs_det(A.cast<lrational_type>().eval()));
            // lint_type det = ModuleNTL::det(A.cast<lint_type>().eval());
            // lint_type det;
            // ModuleFlint::compute_det(A.cast<lint_type>().eval(), det);
            // //  lint_type det = Utility::LinBox::compute_det(A.cast<lint_type>().eval());
            // out() << "determinant computation time: " << timer.finish_str()
            //       << "determinant is: " << det << endl;
        }

        rs_mat solve_AX_eq_B(logger &log, const rs_mat &B)
        {
            if (!decomp_LU)
            {
                decomp_LU = std::make_unique<decomp_type>(A);
            }

            return decomp_LU->solve(B);
        }

        // rs_mat solve_ATX_eq_B(logger &log, const rs_mat &B)
        // {
        //     if (!decomp_LU)
        //     {
        //         decomp_LU = std::make_unique<decomp_type>(A);
        //     }

        //     return decomp_LU->transpose()->solve(B);
        // }

        // size_t get_num_constr() const { return A.rows(); }
    };
}