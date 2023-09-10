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

    namespace Internal
    {
        vector<cut_info> generate_cuts_from_row_matrix__pure_int(const is_mat &cuts_mat, const local_instance &LocalInstance,
                                                                 double tolerance, double err_tolerance)
        {
            RatToInt_Unary rat_to_int(err_tolerance);
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
    }
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

        void append_cuts_from_row_matrices__pure_int(std::vector<is_mat> mat_pool)
        {
            Timer timer;
            out() << "Transforming rays to cuts..." << endl;
            timer.start();

            for (const auto &mat : mat_pool)
            {
                append_cuts(Internal::generate_cuts_from_row_matrix__pure_int(mat.transpose(), inst(),
                                                                              tolerance, err_tolerance));
            }

            out() << "Transforming is finished, time = " << timer.finish_str()
                  << "Filtering duplicate cuts..." << endl;
        }

    public:
        using base::base;

        void generate_cuts() override
        {
            cone_coverer_ptr->direction_strategies.push_back(std::make_unique<typename DirectionStrategies::SpecifiedDirection>(inst().c));

            auto mat_pool = cone_coverer_ptr->generate_cuts(inst().A.transpose().unaryExpr(rat_to_int));

            cone_coverer_ptr->direction_strategies.pop_back();

            append_cuts_from_row_matrices__pure_int(mat_pool);
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
