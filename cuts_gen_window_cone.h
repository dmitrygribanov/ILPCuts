#pragma once

#include "cuts_gen_cone_coverer.h"

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

    class RandomNondegSubmatGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        std::mt19937 &rgen() const { return *_cuts_gen_ref.rgen; }
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        RandomNondegSubmatGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::tuple<std::vector<size_t>, std::vector<size_t>> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<size_t> one_to_n_idxs(n);
            std::iota(one_to_n_idxs.begin(), one_to_n_idxs.end(), 0);

            std::vector<size_t> rows;
            std::sample(one_to_n_idxs.begin(), one_to_n_idxs.end(), std::back_inserter(rows), submat_n, rgen());
            std::sort(rows.begin(), rows.end());

            std::vector<size_t> cols = Utility::generate_nondegenerate_from_rows(inst().A.cast<rational_type>().eval(), rows);

            // size_t rr = Utility::rank(Utility::generate_submatrix_byidx(inst().A.cast<rational_type>().eval(), rows, cols));

            // auto B = Utility::generate_submatrix_byidx(inst.A, rows, cols);

            // if (Utility::rank(B) == submat_n)
            //     return std::make_tuple(rows, cols);

            return std::make_tuple(rows, cols);
        }
    };
    class SomeNondegSubmatGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        SomeNondegSubmatGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::tuple<std::vector<size_t>, std::vector<size_t>> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<size_t> rows(submat_n);
            std::iota(rows.begin(), rows.end(), 0);

            std::vector<size_t> cols = Utility::generate_nondegenerate_from_rows(inst().A, rows);

            return std::make_tuple(rows, cols);
        }
    };

    template <typename SubmatSizeFunT = DimSizeWindow, typename SubmatGenT = RandomNondegSubmatGen>
    class WindowConeCoverCutsGenerator : public ConeCoverCutsGenerator
    {
        using base = ConeCoverCutsGenerator;

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        WindowConeCoverCutsGenerator(logger &log)
            : ConeCoverCutsGenerator(log), submat_gen(*this)
        {
        }

        void generate_cuts() override
        {
            std::vector<size_t> rows, cols;
            std::tie(rows, cols) = submat_gen(submat_size_fun);

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
        }
    };

    using CutsGen_CCRWD = WindowConeCoverCutsGenerator<DimSizeWindow, RandomNondegSubmatGen>;
    using CutsGen_CCRWSqD = WindowConeCoverCutsGenerator<DimSizeWindow, SqrtDimSizeWindow>;
    using CutsGen_CCSWD = WindowConeCoverCutsGenerator<SomeNondegSubmatGen, RandomNondegSubmatGen>;
    using CutsGen_CCSWSqD = WindowConeCoverCutsGenerator<SomeNondegSubmatGen, SqrtDimSizeWindow>;

} // namespace ILPCuts
