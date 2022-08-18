#pragma once

#include "conic_cuts_base.h"
#include "Hilbert_base_alg.h"

// all optional modules are added at the moment
// to simplify experiments
#include "Eigen_solving_utilities.h"
#include "flint_module.h"
#include "ntl_module.h"

#include <omp.h>
#include <atomic>

namespace ILPCuts
{
    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class cuts_generator_extbase
        : public cuts_generator_base<intT, rationalT>
    {
    private:
        using my_type = cuts_generator_extbase<intT, rationalT, lintT, lrationalT>;
        using base = cuts_generator_base<intT, rationalT>;

    public:
        using typename base::int_type;
        using lint_type = lintT;
        using lrational_type = lrationalT;
        using typename base::rational_type;

        using typename base::is_mat;
        using lis_mat = SparseMatrix<lint_type>;
        using lrs_mat = SparseMatrix<lrational_type>;
        using typename base::rs_mat;

        using typename base::is_vec;
        using lis_vec = SparseVector<lint_type>;
        using lrs_vec = SparseVector<lrational_type>;
        using typename base::rs_vec;

        using LRatToLInt_Unary = Utility::RatToInt_Unary<lint_type, lrational_type>;
        LRatToLInt_Unary _lrat_to_lint;
        using base::_rat_to_int;
        static constexpr auto _to_single_int = Utility::conv<lint_type, int_type>;
        static constexpr auto _to_single_rat = Utility::conv<lrational_type, rational_type>;
        using base::_resulting_cuts;
        using base::base;
        using base::cone_limit_per_direction;
        using base::err_tolerance;
        using base::log;
        using base::out;
        using base::rgen;
        using base::threads;
        using base::tolerance;

        using typename base::cut_info;
        using typename base::Strategies;

        using base::copy_original_problem_info;
        using base::generate_cuts_from_row_matrix;
        using base::generate_directions_pool;

        using base::clear_cuts;
        using base::direction_strategies;

    private:
        Hilbert_base_alg<intT, rationalT, lintT, lrationalT> _Hilbert_base_alg;

        void init_Hilbert_alg_params()
        {
            //_Hilbert_base_alg.tolerance = tolerance;
            _Hilbert_base_alg.err_tolerance = err_tolerance;
            //_Hilbert_base_alg.rgen = rgen;
            //_Hilbert_base_alg.threads = threads;
        }

    protected:
        using transform_base = ILPCuts::transform_base<int_type, rational_type, lint_type, lrational_type>;
        using sparse_transform = ILPCuts::sparse_transform<int_type, rational_type, lint_type, lrational_type>;
        using lsparse_transform = ILPCuts::lsparse_transform<int_type, rational_type, lint_type, lrational_type>;
        struct dir_coeffs
        {
            rs_vec coeffs;
            size_t dir_ID;
        };

        struct cone_traits
        {
            is_mat A;
            lint_type det;
            std::vector<dir_coeffs> dir_coeffs_pull;

            void move_from_base_ptr(std::unique_ptr<cone_traits> base)
            {
                //*this = std::move(*base);
                A = std::move(base->A);
                det = std::move(base->det);
                dir_coeffs_pull = std::move(base->dir_coeffs_pull);
            }
        };

        virtual std::unique_ptr<cone_traits> compute_cone_traits(const cone_traits *C,
                                                                 const transform_base &trans,
                                                                 std::vector<dir_coeffs> dir_coeffs_pull)
        {

            size_t n = C->A.cols();
            size_t i = trans.get_i();

            auto new_cone = std::make_unique<cone_traits>();

            new_cone->dir_coeffs_pull = std::move(dir_coeffs_pull);

            try
            {
                out() << "Updating " << i << "-th column of A ...\n";
                new_cone->A = C->A;

                new_cone->A.col(i) = trans(C->A);

                Utility::prune_zeroes_int(new_cone->A);
            }
            catch (ToleranceErrorException exc)
            {
                out() << "A-column construction tolerance error: " << exc.what() << endl;
            }

            try
            {
                out() << "Updating a determinant ...\n";
                new_cone->det = trans[C->det];
            }
            catch (ToleranceErrorException exc)
            {
                out() << "Determinant recompute tolerance error: " << exc.what() << endl;
            }

            return std::move(new_cone);
        }

        std::vector<is_mat> partial_uni_decomp(std::unique_ptr<const cone_traits> C);

        std::vector<std::unique_ptr<cone_traits>> compute_layer(const cone_traits *C,
                                                                std::vector<Atomic<size_t>> &directions_evaluation_counter);

        std::vector<dir_coeffs> compute_dirr_coeffs_pull(const cone_traits &C,
                                                         const transform_base &trans,
                                                         const std::vector<Atomic<size_t>> &directions_counter)
        {

            size_t n = C.A.cols();

            std::vector<dir_coeffs> res;

            for (size_t i = 0; i < C.dir_coeffs_pull.size(); ++i)
            {
                size_t dir_ID = C.dir_coeffs_pull[i].dir_ID;
                if (directions_counter[dir_ID] < cone_limit_per_direction)
                {

                    const rs_vec &old_coeffs = C.dir_coeffs_pull[i].coeffs;
                    rs_vec new_coeffs = trans.compute_coeffs(old_coeffs);

                    if (new_coeffs.size() == n)
                    {
                        res.emplace_back(dir_coeffs{
                            std::move(new_coeffs), dir_ID});
                    }
                }
            }

            return std::move(res);
        }

        using process_func = std::function<void(const transform_base &)>;

        virtual void process_node(const cone_traits &C, process_func pf) = 0;

        void generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP);

    public:
        using size_func = std::function<size_t(size_t)>;
        size_func determinant_cut_limit = [](size_t n) -> size_t
        { return n; };

        explicit cuts_generator_extbase(logger &log)
            : base(log), _Hilbert_base_alg(log), _lrat_to_lint(err_tolerance)
        {
        }
    };

    template <typename intT, typename rationalT, typename lintT, typename lratT>
    std::vector<typename cuts_generator_extbase<intT, rationalT, lintT, lratT>::is_mat>
    cuts_generator_extbase<intT, rationalT, lintT, lratT>::partial_uni_decomp(std::unique_ptr<const cone_traits> C)
    {
        size_t n = C->A.cols();
        omp_set_dynamic(std::max(threads, size_t(1)));
        std::vector<std::unique_ptr<const cone_traits>> current_layer;
        std::vector<std::vector<is_mat>> result_by_treads(threads);

        size_t directions_pool_size = C->dir_coeffs_pull.size();
        current_layer.push_back(std::move(C));

        while (!current_layer.empty())
        {
            std::vector<Atomic<size_t>> directions_evaluation_counter(directions_pool_size);

            size_t dynamic_threads = std::min(threads, current_layer.size());

            std::vector<std::vector<std::unique_ptr<const cone_traits>>> next_layers_by_treads(dynamic_threads);
            omp_set_num_threads(dynamic_threads);
#ifndef ILPCUTS_PARALLEL_OFF
#pragma omp parallel for
#endif
            for (size_t i = 0; i < current_layer.size(); ++i)
            {
                auto node = current_layer[i].get();
                out() << "Computing next layer ..." << endl;
                auto next_node_layer = compute_layer(node, directions_evaluation_counter);
                int tid = omp_get_thread_num();
                for (auto &next_node : next_node_layer)
                {
                    if (next_node->det <= determinant_cut_limit(n))
                    {
                        out() << "Determinant is cutted on the value of " << determinant_cut_limit(n) << ",\n"
                              << "Total Gomory's cuts procedure is started..." << endl;
                        result_by_treads[tid].push_back(_Hilbert_base_alg.compute(next_node->A, next_node->det));
                    }
                    else
                        next_layers_by_treads[tid].push_back(std::move(next_node));
                }
            }

            std::swap(current_layer, next_layers_by_treads[0]);
            for (int i = 1; i < next_layers_by_treads.size(); ++i)
            {
                current_layer.insert(
                    current_layer.end(),
                    std::make_move_iterator(next_layers_by_treads[i].begin()),
                    std::make_move_iterator(next_layers_by_treads[i].end()));
            }
        }

        std::vector<is_mat> &concatenated_result = result_by_treads[0];
        for (int i = 1; i < result_by_treads.size(); ++i)
        {
            concatenated_result.insert(
                concatenated_result.end(),
                std::make_move_iterator(result_by_treads[i].begin()),
                std::make_move_iterator(result_by_treads[i].end()));
        }

        return std::move(concatenated_result);
    }

    template <typename intT, typename rationalT, typename lintT, typename lratT>
    std::vector<std::unique_ptr<typename cuts_generator_extbase<intT, rationalT, lintT, lratT>::cone_traits>>
    cuts_generator_extbase<intT, rationalT, lintT, lratT>::compute_layer(const cone_traits *C,
                                                                         std::vector<Atomic<size_t>> &directions_evaluation_counter)
    {

        std::vector<std::unique_ptr<cone_traits>> result;

        const auto &A = C->A;
        const auto &det = C->det;
        size_t n = A.cols();

        out() << "Current determinant: " << det << "\n"
              << "Nonzeros in A: " << A.nonZeros() << "\n"
              << "A's inifinity norm: " << Utility::max_norm(A) << endl;

        auto process_new_candidate = [&](const transform_base &trans)
        {
            // out() << "Processing new cone candidate ...\n"
            // 	  << "Computing new directions coeff. pull ... \n";
            auto dirr_coeffs_pull = std::move(compute_dirr_coeffs_pull(
                *C, trans, directions_evaluation_counter));

            if (!dirr_coeffs_pull.empty())
            {
                // out() << "New candidate cone is approved ...\n";
                ++directions_evaluation_counter[dirr_coeffs_pull[0].dir_ID];

                std::unique_ptr<cone_traits> new_cone = std::move(compute_cone_traits(
                    C, trans, std::move(dirr_coeffs_pull)));

                result.push_back(std::move(new_cone));
            }
        };

        process_node(*C, process_new_candidate);

        return std::move(result);
    }

    template <typename intT, typename rationalT, typename lintT, typename lratT>
    void cuts_generator_extbase<intT, rationalT, lintT, lratT>::generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP)
    {

        init_Hilbert_alg_params();

        clear_cuts();

        copy_original_problem_info(C->A, optLP);

        Timer timer;
        timer.start();

        out() << "Cone decomposition procedure is started ..." << endl;

        std::vector<is_mat> uni_mat_pool = std::move(partial_uni_decomp(
            std::move(C)));

        out() << "Cone decomposition procedure is finished, time = " << timer.finish_str();

        for (is_mat &uni_mat : uni_mat_pool)
        {
            uni_mat = uni_mat.transpose().eval();
        }

        out() << "Transforming cones to cuts..." << endl;
        timer.start();

        for (const is_mat &uni_mat : uni_mat_pool)
        {
            auto new_cuts = generate_cuts_from_row_matrix(uni_mat, optLP);

            _resulting_cuts.insert(_resulting_cuts.end(), make_move_iterator(new_cuts.begin()),
                                   make_move_iterator(new_cuts.end()));
        }

        out() << "Transforming is finished, time = " << timer.finish_str()
              << "Filtering duplicate cuts..." << endl;

        std::sort(_resulting_cuts.begin(), _resulting_cuts.end());
        auto last = std::unique(_resulting_cuts.begin(), _resulting_cuts.end());
        _resulting_cuts.erase(last, _resulting_cuts.end());
    }
} // namespace ILPCuts
