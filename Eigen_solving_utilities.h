#pragma once

#include "dependencies/Eigen/SparseLU"
#include "dependencies/Eigen/SparseQR"

#include <unordered_map>

#include "basic_utils.h"

namespace ILPCuts
{
    namespace Utility
    {
        class EigenSolveException : public Exception
        {
        public:
            Eigen::ComputationInfo comp_info;

            EigenSolveException(Eigen::ComputationInfo a_comp_info)
                : comp_info(a_comp_info), Exception("ILPCuts lib error: unable to solve/decompose sparse system using Eigen. \n CompInfo: " + std::to_string(a_comp_info)) {}

            EigenSolveException(std::string msg, Eigen::ComputationInfo a_comp_info)
                : comp_info(a_comp_info), Exception(msg) {}
        };

        template <typename ratMatrixT, typename solverT>
        ratMatrixT compute_inverse(const ratMatrixT &A, solverT &solver);

        template <typename rationalT, typename solverT>
        SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A, solverT &solver);

        template <typename rationalT>
        std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A);

        template <typename rationalT>
        rationalT compute_abs_det(const SparseMatrix<rationalT> &A);

        template <typename rationalT>
        SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b);

        template <typename rationalT, typename T2>
        SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, const SparseSolverBase<T2> &solver);

        template <typename intT>
        intT vec_max2pow(const SparseVector<intT> &vec);

        template <>
        mpz_class vec_max2pow(const SparseVector<mpz_class> &vec);

        template <>
        mpz_class vec_max2pow(const SparseVector<mpz_class> &vec)
        {
            using int_type = mpz_class;
            using is_vec = SparseVector<mpz_class>;

            int_type res = 0;

            mp_bitcnt_t min_bit_index = std::numeric_limits<mp_bitcnt_t>::max();
            for (typename is_vec::InnerIterator it(vec); it; ++it)
            {
                if (it.value() != 0)
                {
                    mp_bitcnt_t bit_index = mpz_scan1(it.value().get_mpz_t(), 0);

                    // std::cout << it.index() <<", " << it.value() << std::endl;
                    min_bit_index = std::min(min_bit_index, bit_index);
                }
            }

            mpz_setbit(res.get_mpz_t(), min_bit_index);

            return std::move(res);
        }

        template <typename rationalT, typename solverT>
        SparseMatrix<rationalT> compute_inverse(solverT &solver)
        {
            using rat_mat = SparseMatrix<rationalT>;

            size_t n = solver.cols();

            rat_mat I(n, n);
            I.reserve(VectorXi::Constant(n, 2));
            for (size_t i = 0; i < n; ++i)
                I.insert(i, i) = 1;

            return std::move(solver.solve(I));
        }

        template <typename rationalT>
        SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A)
        {
            SparseLU<SparseMatrix<rationalT>, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            return compute_inverse<rationalT, decltype(decomp_LU)>(decomp_LU);
        }

        template <typename rationalT>
        std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A)
        {

            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            rat_type det = abs(decomp_LU.determinant());
            rat_mat adjA = compute_inverse<rat_type>(A, decomp_LU) * det;

            return std::make_tuple(std::move(adjA), std::move(det));
        }

        template <typename rationalT>
        rationalT compute_abs_det(const SparseMatrix<rationalT> &A)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            rat_type det = abs(decomp_LU.determinant());

            return std::move(det);
        }

        template <typename rationalT>
        size_t rank(const SparseMatrix<rationalT> &A)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseQR<rat_mat, COLAMDOrdering<int>> decomp_QR(A);

            if (decomp_QR.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_QR.info());

            return decomp_QR.rank();
        }

        template <typename rationalT>
        SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b)
        {
            using rat_type = rationalT;
            using rat_mat = SparseMatrix<rat_type>;

            SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

            if (decomp_LU.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_LU.info());

            return std::move(system_solve(b, decomp_LU));
        }

        template <typename rationalT, typename T2>
        SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, const SparseSolverBase<T2> &solver)
        {
            return std::move(solver.solve(b));
        }

        template <typename T>
        SparseMatrix<T> generate_submatrix_byidx(const SparseMatrix<T> &A, const std::vector<size_t> &rows, const std::vector<size_t> &cols)
        {
            using mat_type = SparseMatrix<T>;
            using triplet = Eigen::Triplet<T>;

            size_t m = rows.size();
            size_t n = cols.size();

            assert(n <= A.cols());
            assert(m <= A.rows());

            std::unordered_map<size_t, size_t> rows_hmap(2 * m);
            std::unordered_map<size_t, size_t> cols_hmap(2 * n);

            for (size_t i = 0; i < m; ++i)
            {
                rows_hmap[rows[i]] = i;
            }

            for (size_t j = 0; j < n; ++j)
            {
                cols_hmap[cols[j]] = j;
            }

            mat_type R(m, n);
            std::vector<triplet> tripletsR;

            for (int k = 0; k < A.outerSize(); ++k)
                for (typename mat_type::InnerIterator it(A, k); it; ++it)
                {
                    auto row_hash_pos = rows_hmap.find(it.row());
                    auto col_hash_pos = cols_hmap.find(it.col());
                    if (row_hash_pos != rows_hmap.end() && col_hash_pos != cols_hmap.end())
                        tripletsR.push_back(triplet(row_hash_pos->second, col_hash_pos->second, it.value()));
                }

            R.setFromTriplets(tripletsR.begin(), tripletsR.end());

            // for (size_t i = 0; i < m; ++i)
            // {
            //     for (size_t j = 0; j < n; ++j)
            //     {
            //         R.insert(i, j) = A.coeff(rows[i], cols[j]);
            //     }
            // }

            R.makeCompressed();

            return R;
        }

        template <typename T>
        SparseMatrix<T> generate_row_submatrix_byidx(const SparseMatrix<T> &A, const std::vector<size_t> &rows)
        {
            assert(rows.size() <= A.rows());

            std::vector<size_t> cols(A.cols());
            std::iota(cols.begin(), cols.end(), 0);

            return generate_submatrix_byidx(A, rows, cols);
        }

        template <typename T>
        SparseMatrix<T> generate_column_submatrix_byidx(const SparseMatrix<T> &A, const std::vector<size_t> &cols)
        {
            assert(cols.size() <= A.cols());

            std::vector<size_t> rows(A.rows());
            std::iota(rows.begin(), rows.end(), 0);

            return generate_submatrix_byidx(A, rows, cols);
        }

        template <typename T>
        SparseMatrix<T> generate_subvector_byidx(const SparseVector<T> &v, const std::vector<size_t> &idxs)
        {
            assert(idxs.size() <= v.size());

            SparseVector<T> res(idxs.size());

            for (size_t i = 0; i < idxs.size(); ++i)
            {
                res.insert(i) = v.coeff(idxs[i]);
            }

            return res;
        }

        template <typename T>
        std::vector<size_t> generate_nondegenerate_from_rows(const SparseMatrix<T> &A, std::vector<size_t> rows)
        {
            size_t n = A.cols();
            size_t m = rows.size();

            assert(n == A.rows());
            assert(m <= n);

            using rat_type = T;
            using rs_mat = SparseMatrix<rat_type>;

            rs_mat subA = generate_row_submatrix_byidx(A, rows);

            SparseQR<rs_mat, COLAMDOrdering<int>> decomp_QR(subA);

            if (decomp_QR.info() != ComputationInfo::Success)
                throw EigenSolveException(decomp_QR.info());

            std::vector<size_t> cols(m);

            auto column_perm = decomp_QR.colsPermutation().indices().array();

            // size_t col_perm_size = column_perm.size();
            // size_t a = subA.rows();
            // size_t aa = subA.cols();
            // size_t rr = decomp_QR.rank();

            assert(m <= column_perm.size());

            for (size_t i = 0; i < m; ++i)
            {
                cols[i] = column_perm[i];
            }

            std::sort(cols.begin(), cols.end());

            // auto subsubA = generate_submatrix_byidx(A, rows, cols);

            // decomp_QR.compute(subsubA);
            // rr = decomp_QR.rank();
            // size_t rrr = rank(subsubA);

            return cols;
        }
    } // namespace Utility

} // namespace ILPCuts
