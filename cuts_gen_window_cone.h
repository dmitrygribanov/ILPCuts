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

    class RandomNondegSubmatGen
    {
        cuts_gen_base &_cuts_gen_ref;

    public:
        std::mt19937 &rgen() const { return *_cuts_gen_ref.rgen; }
        const local_instance &inst() const { return _cuts_gen_ref.LocalInstance; }

        RandomNondegSubmatGen(cuts_gen_base &cuts_gen_ref)
            : _cuts_gen_ref(cuts_gen_ref) {}

        template <typename submatSizeFunT>
        std::vector<size_t> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<size_t> one_to_n_idxs(n);
            std::iota(one_to_n_idxs.begin(), one_to_n_idxs.end(), 0);

            std::vector<size_t> rows;
            std::sample(one_to_n_idxs.begin(), one_to_n_idxs.end(), std::back_inserter(rows), submat_n, rgen());
            std::sort(rows.begin(), rows.end());

            return rows;
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
        std::vector<size_t> operator()(submatSizeFunT submat_size_fun) const
        {
            size_t n = inst().get_num_vars();
            size_t submat_n = submat_size_fun(inst());

            std::vector<size_t> rows(submat_n);
            std::iota(rows.begin(), rows.end(), 0);

            return rows;
        }
    };

    template <typename SubmatSizeFunT = DimSizeWindow, typename SubmatGenT = RandomNondegSubmatGen>
    class ___WindowConeCoverCutsGenerator : public ConeCoverCutsGenerator
    {
        using base = ConeCoverCutsGenerator;

    public:
        SubmatSizeFunT submat_size_fun;
        SubmatGenT submat_gen;

        ___WindowConeCoverCutsGenerator(logger &log)
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
            Timer timer;
            out() << "Window cutting plane generation is started..." << endl;

            std::vector<size_t> rows;
            rows = submat_gen(submat_size_fun);

            is_mat rowSubA = Utility::generate_row_submatrix_byidx(LocalInstance.A, rows);
            out() << "Row submatrix generation is finished..." << endl;

            is_mat Q, H;

            timer.start();
            out() << "HNF computation is started..." << endl;
            ModuleFlint::compute_HNF<int_type>(rowSubA.transpose(), Q, H);
            Q = Q.transpose().eval();
            H = H.transpose().eval();
            out() << "HNF computation is finished,\n"
                  << "time: " << timer.finish_str() << endl;

            // is_mat Qi = Utility::compute_inverse(Q.cast<rational_type>().eval()).unaryExpr(rat_to_int);

            std::vector<size_t> cols = Utility::get_indices_of_nonzero_columns(H);

            size_t n = H.cols();
            size_t m = H.rows();
            size_t cc = cols.size();

            is_mat subH = Utility::generate_column_submatrix_byidx(H, cols);

            timer.start();
            out() << "c-vector projecting is satrted..." << endl;
            rs_vec t = Utility::system_solve<rational_type>(LocalInstance.A.cast<rational_type>().transpose(), LocalInstance.c);
            rs_vec t_proj = Utility::generate_subvector_byidx(t, rows);
            rs_vec c_proj = rowSubA.cast<rational_type>().transpose() * t_proj;

            rs_vec c_trans = Q.cast<rational_type>().transpose() * c_proj;
            c_trans = Utility::generate_subvector_byidx(c_trans, cols);
            out() << "c-vector projecting is finished,\n"
                  << "time: " << timer.finish_str() << endl;

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

            out() << "Window cutting plane generation is finished..." << endl;
        }
    };

    using CutsGen_CCRWD = WindowConeCoverCutsGenerator<DimSizeWindow, RandomNondegSubmatGen>;
    using CutsGen_CCRWSqD = WindowConeCoverCutsGenerator<DimSizeWindow, SqrtDimSizeWindow>;
    using CutsGen_CCSWD = WindowConeCoverCutsGenerator<SomeNondegSubmatGen, RandomNondegSubmatGen>;
    using CutsGen_CCSWSqD = WindowConeCoverCutsGenerator<SomeNondegSubmatGen, SqrtDimSizeWindow>;

} // namespace ILPCuts
