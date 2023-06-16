#pragma once

#include "cuts_gen_base.h"
#include "cone_coverer_base.h"

namespace ILPCuts
{
    class ConeCoverCutsGenerator : public cuts_gen_base
    {
    private:
        using my_type = ConeCoverCutsGenerator;
        using base = cuts_gen_base;

    public:
        std::unique_ptr<ConeCoverer_Base> cone_coverer_ptr;

    protected:
        vector<cut_info> generate_cuts_from_row_matrix(const is_mat &cuts_mat)
        {
            vector<cut_info> res;

            rs_vec rat_b = cuts_mat.template cast<rational_type>() * LocalInstance.v;
            is_vec rounded_b = rat_b.unaryExpr(Utility::Round_Unary<rational_type>()).unaryExpr(_rat_to_int);
            is_vec floor_b = rat_b.unaryExpr(Utility::Floor_Unary<rational_type>()).unaryExpr(_rat_to_int);

            rs_vec diff_b = (rat_b - rounded_b.template cast<rational_type>()).cwiseAbs();

            typename rs_vec::InnerIterator diff_it(diff_b);
            while (diff_it)
            {
                if (diff_it.value() > rational_type(tolerance))
                {
                    is_vec ray = cuts_mat.row(diff_it.index());
                    int_type rhs = floor_b.coeff(diff_it.index());
                    rational_type frac_rhs = rat_b.coeff(diff_it.index());

                    res.emplace_back(cut_info{ray, rhs, frac_rhs});
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

            for (auto &mat : mat_pool)
            {
                mat = mat.transpose().eval();
            }

            for (const auto &mat : mat_pool)
            {
                auto new_cuts = generate_cuts_from_row_matrix(mat);

                _resulting_cuts.insert(_resulting_cuts.end(), make_move_iterator(new_cuts.begin()),
                                       make_move_iterator(new_cuts.end()));
            }

            out() << "Transforming is finished, time = " << timer.finish_str()
                  << "Filtering duplicate cuts..." << endl;

            std::sort(_resulting_cuts.begin(), _resulting_cuts.end());
            auto last = std::unique(_resulting_cuts.begin(), _resulting_cuts.end());
            _resulting_cuts.erase(last, _resulting_cuts.end());
        }

    public:
        using base::base;

        void generate_cuts() override
        {
            auto mat_pool = cone_coverer_ptr->generate_cuts(LocalInstance.A.transpose());
            append_cuts_from_row_matrices(mat_pool);
        }

        void set_params(double a_tolerance, double a_err_tolerance, size_t a_threads, std::mt19937 *a_rgen)
        {
            cone_coverer_ptr->tolerance = tolerance = a_tolerance;
            cone_coverer_ptr->err_tolerance = err_tolerance = a_err_tolerance;
            cone_coverer_ptr->threads = threads = a_threads;
            cone_coverer_ptr->rgen = rgen = a_rgen;
        }
    };
} // namespace ILPCuts
