#pragma once

#include "cuts_gen_cone_coverer.h"
#include "flint_module.h"

namespace ILPCuts
{
    class ConstantSizeGen
    {
    public:
        size_t wSize = 0;

        size_t operator()(const local_instance &) const { return wSize; }
    };

    class DimSizeGen
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

    class SqrtDimSizeGen
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

    class RandomIntVarsColsGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        size_t trials_num = 1;

        std::mt19937 &rgen() const { return *_cuts_gen_ref.rgen; }
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        RandomIntVarsColsGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<std::vector<size_t>> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<std::vector<size_t>> res;

            std::vector<size_t> int_vars_idxs_pool;
            for (typename rs_vec::InnerIterator it(inst().v); it; ++it)
            {
                if (is_fract(it.value(), _cuts_gen_ref.tolerance) &&
                    inst().ILP_var_types[it.index()] == local_instance::var_ILP_type::Int)
                {
                    int_vars_idxs_pool.push_back(it.index());
                }
            }
            if (int_vars_idxs_pool.size() < submat_n)
            {
                for (typename rs_vec::InnerIterator it(inst().v); it && int_vars_idxs_pool.size() < submat_n; ++it)
                {
                    if (!is_fract(it.value(), _cuts_gen_ref.tolerance) &&
                        inst().ILP_var_types[it.index()] == local_instance::var_ILP_type::Int)
                    {
                        int_vars_idxs_pool.push_back(it.index());
                    }
                }
            }

            for (size_t i = 0; i < trials_num; ++i)
            {
                std::vector<size_t> cols;
                std::sample(int_vars_idxs_pool.begin(), int_vars_idxs_pool.end(), std::back_inserter(cols), submat_n, rgen());
                std::sort(cols.begin(), cols.end());

                res.push_back(cols);
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

                std::vector<size_t> cols = Utility::generate_nondegenerate_from_rows(inst().A, rows);

                res.push_back(std::make_tuple(rows, cols));
            }

            return res;
        }
    };

    class RandomIntVarsSubmatGen : public RandomIntVarsColsGen
    {
    public:
        RandomIntVarsSubmatGen(cuts_gen_base &cuts_gen_ref)
            : RandomIntVarsColsGen(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<submat_idxs> operator()(submatSizeFunT submat_size_fun) const
        {
            auto cols_trials = RandomIntVarsColsGen::operator()<submatSizeFunT>(submat_size_fun);

            std::vector<submat_idxs> res;

            for (size_t i = 0; i < cols_trials.size(); ++i)
            {
                const std::vector<size_t> &cols = cols_trials[i];

                std::vector<size_t> rows = Utility::generate_nondegenerate_from_rows<rational_type>(inst().A.transpose(), cols);

                res.push_back(std::make_tuple(rows, cols));
            }

            return res;
        }
    };

    namespace Internal
    {
        std::vector<cut_info> generate_split_cuts_from_base_ineq(const local_instance &inst, const is_mat &base_ineq_cols,
                                                                 double tolerance, double err_tolerance)
        {
            std::vector<cut_info> res_cuts;
            RatToInt_Unary rat_to_int(err_tolerance);

            size_t n = inst.get_num_vars();
            assert(n == base_ineq_cols.rows());
            // rs_vec part_v = Utility::generate_subvector_byidx(inst.v, intComp_idxs);

            // assert(part_v.size() == base_ineq_cols.rows());

            rs_vec rat_b = base_ineq_cols.cast<rational_type>().transpose() * inst.v;
            is_vec rounded_b = rat_b.unaryExpr(Utility::Round_Unary<rational_type>()).unaryExpr(rat_to_int);
            is_vec floor_b = rat_b.unaryExpr(Utility::Floor_Unary<rational_type>()).unaryExpr(rat_to_int);
            // rs_vec ceil_b = rat_b.unaryExpr(Utility::Ceil_Unary<rational_type>());

            rs_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).cwiseAbs();
            rs_vec fract_part_b = rat_b - floor_b.cast<rational_type>();

            SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU(-inst.A.transpose());
            if (decomp_LU.info() != ComputationInfo::Success)
                throw Utility::EigenSolveException(decomp_LU.info());

            rs_mat cone_ray_intersections = decomp_LU.solve(base_ineq_cols.cast<rational_type>());

            typename rs_vec::InnerIterator diff_it(diff_b);
            while (diff_it)
            {
                if (!is_zero(diff_it.value(), tolerance))
                {
                    rational_type epsilon = fract_part_b.coeff(diff_it.index());
                    rs_vec intersections = cone_ray_intersections.col(diff_it.index());

                    rs_vec new_ray(n);
                    rational_type new_rhs = 0.0;

                    for (typename rs_vec::InnerIterator it(intersections); it; ++it)
                    {
                        if (!is_zero(it.value(), tolerance))
                        {
                            rational_type alpha = 0.0;
                            if (it.value() > 0)
                            {
                                alpha = (1.0 - epsilon) / it.value();
                            }
                            else
                            {
                                // it.value() < 0
                                alpha = (-epsilon) / it.value();
                            }
                            new_ray += inst.A.row(it.index()).transpose() / alpha;
                        }
                    }

                    new_rhs = -1.0 + new_ray.dot(inst.v);

                    std::vector<size_t> int_idxs(n);
                    std::iota(int_idxs.begin(), int_idxs.end(), 0);

                    res_cuts.emplace_back(cut_info{new_ray, new_rhs, new_rhs, int_idxs});
                }
                ++diff_it;
            }

            return res_cuts;
        }
    }

    template <typename SubmatSizeFunT = DimSizeGen, typename SubmatGenT = RandomIntVarsSubmatGen>
    class Split_SubmatCutsGen : public ConeCoverCutsGenerator
    {
    private:
        using base = ConeCoverCutsGenerator;

        is_mat complete_by_zeros(const is_mat &base_mat, const std::vector<size_t> &var_idxs)
        {
            size_t n = LocalInstance.get_num_vars();

            is_mat res(n, base_mat.cols());

            using triplet = Eigen::Triplet<int_type>;

            std::vector<triplet> triplets;

            for (int k = 0; k < base_mat.outerSize(); ++k)
                for (typename is_mat::InnerIterator it(base_mat, k); it; ++it)
                {
                    triplets.push_back(triplet(var_idxs[it.row()], it.col(), it.value()));
                }

            res.setFromTriplets(triplets.begin(), triplets.end());
            res.makeCompressed();

            return res;
        }

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        Split_SubmatCutsGen(logger &log)
            : ConeCoverCutsGenerator(log), submat_gen(*this)
        {
        }

        void generate_cuts() override
        {
            size_t n = LocalInstance.get_num_vars();
            std::vector<submat_idxs> submat_trials = submat_gen(submat_size_fun);

            for (size_t trial_i = 0; trial_i < submat_trials.size(); ++trial_i)
            {
                std::vector<size_t> rows, cols;
                std::tie(rows, cols) = submat_trials[trial_i];

                rs_mat subA = Utility::generate_submatrix_byidx(LocalInstance.A, rows, cols);
                // is_mat rowSubA = Utility::generate_row_submatrix_byidx(LocalInstance.A, rows);

                rs_vec t = Utility::system_solve<rational_type>(LocalInstance.A.transpose(), LocalInstance.c);
                rs_vec t_proj = Utility::generate_subvector_byidx(t, rows);
                rs_vec c_proj = subA.cast<rational_type>().transpose() * t_proj;

                cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename DirectionStrategies::SpecifiedDirection>(c_proj));

                auto subM_pool = cone_coverer_ptr->generate_cuts(subA.transpose().unaryExpr(rat_to_int));

                cone_coverer_ptr->direction_strategies.pop_back();

                //--------------------------------------------------------------------

                SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU_subA(subA.cast<rational_type>().transpose());

                if (decomp_LU_subA.info() != ComputationInfo::Success)
                    throw Utility::EigenSolveException(decomp_LU_subA.info());

                std::vector<rs_mat> res_cuts_pool;
                for (size_t i = 0; i < subM_pool.size(); ++i)
                {
                    // const is_mat &C = subM_pool[i];
                    // rs_mat T = decomp_LU_subA.solve(C.cast<rational_type>()).transpose();

                    // is_mat cutsMat = T * subA.cast<rational_type>();

                    is_mat cutsBaseMat = complete_by_zeros(subM_pool[i], cols);

                    append_cuts(Internal::generate_split_cuts_from_base_ineq(LocalInstance, cutsBaseMat,
                                                                             tolerance, err_tolerance));

                    // res_cuts_pool.push_back(cutsMat_rat.transpose());
                }
            }

            make_cuts_unique();
        }
    };

    template <typename SubmatSizeFunT = DimSizeGen, typename SubmatGenT = RandomNondegRowsGen>
    class Hermite_SubmatCutsGen : public ConeCoverCutsGenerator
    {
        using base = ConeCoverCutsGenerator;

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        Hermite_SubmatCutsGen(logger &log)
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

                is_mat rowSubA = Utility::generate_row_submatrix_byidx(LocalInstance.A, rows).unaryExpr(rat_to_int);
                out() << "Row submatrix generation is finished..." << endl;

                timer.start();
                out() << "c-vector projecting is satrted..." << endl;
                rs_vec t = Utility::system_solve<rational_type>(LocalInstance.A.transpose(), LocalInstance.c);
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

                append_cuts_from_row_matrices__pure_int(res_cuts_pool);

                out() << "Cutting planes extraction is finished,\n"
                      << "time: " << timer.finish_str() << endl;
            }

            make_cuts_unique();
            out() << "Window cutting plane generation is finished..." << endl;
        }
    };

    using HermiteCutsGen_RWDS = Hermite_SubmatCutsGen<DimSizeGen, RandomNondegRowsGen>;
    using HermiteCutsGen_RWSqS = Hermite_SubmatCutsGen<SqrtDimSizeGen, RandomNondegRowsGen>;

    using SplitCutsGen_RWDS = Split_SubmatCutsGen<DimSizeGen, RandomIntVarsSubmatGen>;

} // namespace ILPCuts
