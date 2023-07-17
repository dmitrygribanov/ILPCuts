#pragma once

#include "cuts_gen_base.h"
#include "cone_coverer_base.h"

namespace ILPCuts
{
    // class CutFunctionBase : public base_alg
    // {
    // public:
    //     CutFunctionBase(logger &log)
    //         : base_alg(log) {}

    //     virtual vector<cut_info> operator()(const rs_mat &cuts_mat) = 0;
    // };

    // class GomoryPureInt_CutF : public

    class ConeCoverCutsGenerator : public cuts_gen_base
    {
    private:
        using my_type = ConeCoverCutsGenerator;
        using base = cuts_gen_base;

    public:
        std::unique_ptr<ConeCoverer_Base> cone_coverer_ptr;

    protected:
        void make_cuts_unique()
        {
            std::sort(_resulting_cuts.begin(), _resulting_cuts.end());
            auto last = std::unique(_resulting_cuts.begin(), _resulting_cuts.end());
            _resulting_cuts.erase(last, _resulting_cuts.end());
        }

        void append_cuts(vector<cut_info> new_cuts)
        {
            _resulting_cuts.insert(_resulting_cuts.end(), make_move_iterator(new_cuts.begin()),
                                   make_move_iterator(new_cuts.end()));
        }

        vector<cut_info> generate_cuts_from_row_matrix(const rs_mat &cuts_mat, const std::vector<size_t> &intComp_idxs)
        {
            vector<cut_info> res;

            size_t n = cuts_mat.cols();
            size_t cuts_num = cuts_mat.rows();
            size_t int_n = intComp_idxs.size();
            size_t rat_n = n - int_n;

            std::vector<size_t> ratComp_idxs = Utility::idxs_set_complement(n, intComp_idxs);

            is_mat intCompsMat = Utility::generate_column_submatrix_byidx(cuts_mat, intComp_idxs).unaryExpr(rat_to_int);
            rs_mat ratCompsMat = Utility::generate_column_submatrix_byidx(cuts_mat, ratComp_idxs);
            rs_vec intComp_v = Utility::generate_subvector_byidx(LocalInstance.v, intComp_idxs);
            rs_vec ratComp_v = Utility::generate_subvector_byidx(LocalInstance.v, ratComp_idxs);

            rs_vec b = cuts_mat * LocalInstance.v;
            rs_vec intComp_b = intCompsMat.cast<rational_type>() * intComp_v;
            rs_vec ratComp_b = b - intComp_b;
            rs_vec offset_b = b - ratComp_b;

            return res;
        }

        vector<cut_info> generate_cuts_from_row_matrix(const is_mat &cuts_mat)
        {
            vector<cut_info> res;

            rs_vec rat_b = cuts_mat.cast<rational_type>() * LocalInstance.v;
            is_vec rounded_b = rat_b.unaryExpr(Utility::Round_Unary<rational_type>()).unaryExpr(rat_to_int);
            is_vec floor_b = rat_b.unaryExpr(Utility::Floor_Unary<rational_type>()).unaryExpr(rat_to_int);

            rs_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).cwiseAbs();

            typename rs_vec::InnerIterator diff_it(diff_b);
            while (diff_it)
            {
                if (!is_zero(diff_it.value(), tolerance))
                {
                    is_vec ray = cuts_mat.row(diff_it.index());
                    int_type rhs = floor_b.coeff(diff_it.index());
                    rational_type precut_rhs = rat_b.coeff(diff_it.index());

                    std::vector<size_t> int_idxs(ray.size());
                    std::iota(int_idxs.begin(), int_idxs.end(), 0);

                    res.emplace_back(cut_info{ray.cast<rational_type>(), static_cast<rational_type>(rhs), precut_rhs, int_idxs});
                }
                ++diff_it;
            }

            return res;
        }

        void append_cuts_from_row_matrices(std::vector<is_mat> mat_pool)
        {
            Timer timer;
            out() << "Transforming rays to cuts..." << endl;
            timer.start();

            for (const auto &mat : mat_pool)
            {
                append_cuts(generate_cuts_from_row_matrix(mat.transpose()));
            }

            out() << "Transforming is finished, time = " << timer.finish_str()
                  << "Filtering duplicate cuts..." << endl;
        }

    public:
        using base::base;

        void generate_cuts() override
        {
            cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename DirectionStrategies::SpecifiedDirection>(LocalInstance.c));

            auto mat_pool = cone_coverer_ptr->generate_cuts(LocalInstance.A.transpose());

            cone_coverer_ptr->direction_strategies.pop_back();

            append_cuts_from_row_matrices(mat_pool);
            make_cuts_unique();
        }

        void set_params(double a_tolerance, double a_err_tolerance, size_t a_threads, std::mt19937 *a_rgen)
        {
            cone_coverer_ptr->tolerance = tolerance = a_tolerance;
            cone_coverer_ptr->err_tolerance = err_tolerance = a_err_tolerance;
            cone_coverer_ptr->threads = threads = a_threads;
            cone_coverer_ptr->rgen = rgen = a_rgen;
        }
    };

    using CutsGen_CC = ConeCoverCutsGenerator;
} // namespace ILPCuts
