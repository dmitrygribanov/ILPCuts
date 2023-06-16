#pragma once

#include "base_algo.h"

#include "Eigen_solving_utilities.h"
#include "flint_module.h"

namespace ILPCuts
{
    class Hilbert_base_alg
        : public base_alg
    {
    private:
        using my_type = Hilbert_base_alg;
        using base = base_alg;

    private:
        void fill_parall_vectors(const std::vector<std::vector<int_type>> &A, size_t i, std::vector<rational_type> &x, std::vector<int_type> &y, std::vector<std::vector<int_type>> &R)
        {
            size_t n = x.size();
            rational_type r_i = 0;
            for (size_t k = i + 1; k < n; ++k)
                r_i += x[k] * A[i][k];
            int_type ceil_r_i = conv<rational_type, int_type>::make(ceil(r_i));
            int_type initial_y = ceil_r_i;
            int_type H_ii = A[i][i];
            for (int_type y_i = initial_y; y_i < initial_y + H_ii; ++y_i)
            {
                rational_type x_i = (y_i - r_i) / H_ii;
                x[i] = x_i;
                y[i] = y_i;
                if (i != 0)
                    fill_parall_vectors(A, i - 1, x, y, R);
                else
                    R.emplace_back(y);
            }
        }

        is_mat compute_parall_vectors(const std::vector<std::vector<int_type>> &A)
        {
            size_t n = A.size();
            std::vector<std::vector<int_type>> R_dense;

            std::vector<rational_type> x(n, 0.0);
            std::vector<int_type> y(n, 0);
            fill_parall_vectors(A, n - 1, x, y, R_dense);

            size_t m = R_dense.size();
            is_mat R_sparse(n, m);

            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i = 0; i < n; ++i)
                {
                    if (R_dense[j][i] != 0)
                        R_sparse.insert(i, j) = R_dense[j][i];
                }
            }

            return R_sparse;
        }

    public:
        using base::base;

        logger::stream &out()
        {
            log.out() << "ILPCuts-Hilbert-base-alg: ";
            return log.out();
        }

        is_mat compute(const is_mat &A, const lint_type &det)
        {
            Timer timer;
            size_t n = A.cols();
            is_mat Pi, H;
            rs_mat Hi;
            // int_type denum;

            out() << "computation of the HNF is started...\n";
            timer.start();
            // ModuleFlint::compute_HNF(C.A, C.A);
            // H = ModuleNTL::compute_HNF<int_type>(C.A.transpose(), to_single_int(C.det));
            ModuleFlint::compute_HNF(A, to_single_int(det), H);

            // cout << "Is Upper triang: " << Utility::is_upper_triang(H) << endl;

            out() << "HNF computation time: " << timer.finish_str();

            out() << "computation of the H^-1 is started...\n";
            timer.start();
            Hi = Utility::compute_inverse(H.cast<rational_type>().eval());
            //  ModuleFlint::compute_inverse<int64_t>(H, denum, Hi);
            out() << "HNF^-1 computation time: " << timer.finish_str();

            out() << "computation of the P^-1 is started...\n";
            timer.start();
            Pi = (A.cast<rational_type>() * Hi).unaryExpr(rat_to_int);
            out() << "P^-1 computation time: " << timer.finish_str();

            std::vector<std::vector<int_type>> H_dense;
            for (size_t i = 0; i < n; ++i)
            {
                H_dense.emplace_back(std::vector<int_type>(n, 0));
            }

            for (size_t k = 0; k < H.outerSize(); ++k)
                for (typename is_mat::InnerIterator it(H, k); it; ++it)
                {
                    if (it.value() != 0)
                        H_dense[it.row()][it.col()] = it.value();
                }

            out() << "parall vectors computation is started...\n";
            timer.start();
            is_mat parall_vectors = compute_parall_vectors(H_dense);
            out() << "parall vectors computation time: " << timer.finish_str();

            is_mat res(n, n + parall_vectors.cols());
            res.leftCols(n) = A;
            res.middleCols(n, parall_vectors.cols()) = Pi * parall_vectors;

            return res;
        }
    };
} // namespace ILPCuts
