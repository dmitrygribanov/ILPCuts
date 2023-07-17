#pragma once

#include "cuts_gen_cone_coverer.h"
#include "flint_module.h"

namespace ILPCuts
{
    class ConstantSizeWindow
    {
    public:
        size_t wSize = 0;

        size_t operator()(const local_instance &) const { return wSize; }
    };

    class DimSizeWindow
    {
    public:
        size_t maxSize = 0;

        size_t operator()(const local_instance &inst) const
        {
            size_t n = inst.get_num_vars();

            if (n > maxSize)
                return maxSize;

            return n;
        }
    };

    class SqrtDimSizeWindow
    {
    public:
        size_t maxSize = 0;
        size_t minSize = 0;

        size_t operator()(const local_instance &inst) const
        {
            size_t n = inst.get_num_vars();
            double sqrt_n = std::floor(
                std::sqrt(static_cast<double>(n)));

            if (sqrt_n < minSize)
                return minSize;

            if (sqrt_n > maxSize)
                return maxSize;

            return sqrt_n;
        }
    };

    class RandomNondegRowsGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        size_t trials_num = 1;

        std::mt19937 &rgen() const { return *_cuts_gen_ref.rgen; }
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        RandomNondegRowsGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<std::vector<size_t>> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<std::vector<size_t>> res;

            std::vector<size_t> one_to_n_idxs(n);
            std::iota(one_to_n_idxs.begin(), one_to_n_idxs.end(), 0);

            for (size_t i = 0; i < trials_num; ++i)
            {
                std::vector<size_t> rows;
                std::sample(one_to_n_idxs.begin(), one_to_n_idxs.end(), std::back_inserter(rows), submat_n, rgen());
                std::sort(rows.begin(), rows.end());

                res.push_back(rows);
            }

            return res;
        }
    };

    using submat_idxs = std::tuple<std::vector<size_t>, std::vector<size_t>>;

    class RandomNondegSubmatGen : public RandomNondegRowsGen
    {
    public:
        RandomNondegSubmatGen(cuts_gen_base &cuts_gen_ref)
            : RandomNondegRowsGen(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<submat_idxs> operator()(submatSizeFunT submat_size_fun) const
        {
            auto rows_trials = RandomNondegRowsGen::operator()<submatSizeFunT>(submat_size_fun);

            std::vector<submat_idxs> res;

            for (size_t i = 0; i < rows_trials.size(); ++i)
            {
                const std::vector<size_t> &rows = rows_trials[i];
                is_mat rowSubA = Utility::generate_row_submatrix_byidx(inst().A, rows);

                std::vector<size_t> cols = Utility::generate_nondegenerate_from_rows(inst().A, rows);

                res.push_back(std::make_tuple(rows, cols));
            }

            return res;
        }
    };
    class SomeNondegRowsGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        SomeNondegRowsGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<std::vector<size_t>> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<size_t> rows(submat_n);
            std::iota(rows.begin(), rows.end(), 0);

            std::vector<std::vector<size_t>> res;
            res.push_back(rows);

            return res;
        }
    };

    template <typename SubmatSizeFunT = DimSizeWindow, typename SubmatGenT = RandomNondegSubmatGen>
    class SqSubmat_WindowCCCutsGen : public ConeCoverCutsGenerator
    {
        using base = ConeCoverCutsGenerator;

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        SqSubmat_WindowCCCutsGen(logger &log)
            : ConeCoverCutsGenerator(log), submat_gen(*this)
        {
        }

        void generate_cuts() override
        {
            std::vector<submat_idxs> submat_trials = submat_gen(submat_size_fun);

            for (size_t trial_i = 0; trial_i < submat_trials.size(); ++trial_i)
            {
                std::vector<size_t> rows, cols;
                std::tie(rows, cols) = submat_trials[trial_i];

                is_mat subA = Utility::generate_submatrix_byidx(LocalInstance.A, rows, cols);
                is_mat rowSubA = Utility::generate_row_submatrix_byidx(LocalInstance.A, rows);

                rs_vec t = Utility::system_solve<rational_type>(LocalInstance.A.cast<rational_type>().transpose(), LocalInstance.c);
                rs_vec t_proj = Utility::generate_subvector_byidx(t, rows);
                rs_vec c_proj = subA.cast<rational_type>().transpose() * t_proj;

                cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename DirectionStrategies::SpecifiedDirection>(c_proj));

                auto subM_pool = cone_coverer_ptr->generate_cuts(subA.transpose());

                cone_coverer_ptr->direction_strategies.pop_back();

                //--------------------------------------------------------------------

                SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU_subA(subA.cast<rational_type>().transpose());

                if (decomp_LU_subA.info() != ComputationInfo::Success)
                    throw Utility::EigenSolveException(decomp_LU_subA.info());

                std::vector<rs_mat> res_cuts_pool;
                for (size_t i = 0; i < subM_pool.size(); ++i)
                {
                    const is_mat &C = subM_pool[i];
                    rs_mat T = decomp_LU_subA.solve(C.cast<rational_type>()).transpose();

                    rs_mat cutsMat_rat = T * rowSubA.cast<rational_type>();

                    res_cuts_pool.push_back(cutsMat_rat.transpose());
                }

                // append_cuts_from_row_matrices(res_cuts_pool);
            }

            make_cuts_unique();
        }
    };

    template <typename SubmatSizeFunT = DimSizeWindow, typename SubmatGenT = RandomNondegRowsGen>
    class Hermite_WindowCCCutsGen : public ConeCoverCutsGenerator
    {
        using base = ConeCoverCutsGenerator;

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        Hermite_WindowCCCutsGen(logger &log)
            : ConeCoverCutsGenerator(log), submat_gen(*this)
        {
        }

        void generate_cuts() override
        {
            std::vector<std::vector<size_t>> rows_trials;
            rows_trials = submat_gen(submat_size_fun);

            Timer timer;
            out() << "Window cutting plane generation is started, number of trials: " << rows_trials.size() << endl;

            for (size_t trial_i = 0; trial_i < rows_trials.size(); ++trial_i)
            {
                const auto &rows = rows_trials[trial_i];
                out() << "Trial index: " << trial_i << endl;

                is_mat rowSubA = Utility::generate_row_submatrix_byidx(LocalInstance.A, rows);
                out() << "Row submatrix generation is finished..." << endl;

                timer.start();
                out() << "c-vector projecting is satrted..." << endl;
                rs_vec t = Utility::system_solve<rational_type>(LocalInstance.A.cast<rational_type>().transpose(), LocalInstance.c);
                is_vec t_proj = Utility::generate_subvector_byidx(t, rows).unaryExpr(Utility::Ceil_Unary<rational_type>()).cast<int_type>();
                // rs_vec c_proj = rowSubA.cast<rational_type>().transpose() * t_proj;
                is_vec c_proj = rowSubA.transpose() * t_proj;
                out() << "c-vector projecting is finished,\n"
                      << "time: " << timer.finish_str() << endl;

                timer.start();
                out() << "HNF computation is started..." << endl;
                is_mat H;
                is_vec c_trans;
                ModuleFlint::compute_HNF<int_type>(rowSubA.transpose(), c_proj, H, c_trans);
                H = H.transpose().eval();

                std::vector<size_t> cols = Utility::get_indices_of_nonzero_columns(H);

                is_mat subH = Utility::generate_column_submatrix_byidx(H, cols);
                c_trans = Utility::generate_subvector_byidx(c_trans, cols);

                // rs_vec c_trans = Q.cast<rational_type>().transpose() * c_proj;
                // c_trans = Utility::generate_subvector_byidx(c_trans, cols);

                out() << "HNF computation is finished,\n"
                      << "time: " << timer.finish_str() << endl;

                // is_mat Qi = Utility::compute_inverse(Q.cast<rational_type>().eval()).unaryExpr(rat_to_int);

                cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename DirectionStrategies::SpecifiedDirection>(c_trans));

                auto subM_pool = cone_coverer_ptr->generate_cuts(subH.transpose());

                cone_coverer_ptr->direction_strategies.pop_back();

                //--------------------------------------------------------------------

                timer.start();
                out() << "Cutting planes extraction is started..." << endl;
                SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU_subA(subH.cast<rational_type>().transpose());

                if (decomp_LU_subA.info() != ComputationInfo::Success)
                    throw Utility::EigenSolveException(decomp_LU_subA.info());

                std::vector<is_mat> res_cuts_pool;
                for (size_t i = 0; i < subM_pool.size(); ++i)
                {
                    const is_mat &C = subM_pool[i];
                    rs_mat T = decomp_LU_subA.solve(C.cast<rational_type>()).transpose();

                    rs_mat cutsMat_rat = T * rowSubA.cast<rational_type>();
                    is_mat cutsMat_int = cutsMat_rat.unaryExpr(rat_to_int);

                    res_cuts_pool.push_back(cutsMat_int.transpose());
                }

                append_cuts_from_row_matrices(res_cuts_pool);

                out() << "Cutting planes extraction is finished,\n"
                      << "time: " << timer.finish_str() << endl;
            }

            make_cuts_unique();
            out() << "Window cutting plane generation is finished..." << endl;
        }
    };

    using HermiteCutsGen_RWDS = Hermite_WindowCCCutsGen<DimSizeWindow, RandomNondegRowsGen>;
    using HermiteCutsGen_RWSqS = Hermite_WindowCCCutsGen<SqrtDimSizeWindow, RandomNondegRowsGen>;
    using HermiteCutsGen_SWDS = Hermite_WindowCCCutsGen<DimSizeWindow, SomeNondegRowsGen>;
    using HermiteCutsGen_SWSqS = Hermite_WindowCCCutsGen<SqrtDimSizeWindow, SomeNondegRowsGen>;

} // namespace ILPCuts
