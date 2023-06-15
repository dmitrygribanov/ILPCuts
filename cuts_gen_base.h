#pragma once

#include "cone_coverer_base.h"
#include "local_instance.h"
#include "basic_utils.h"

namespace ILPCuts
{
    template <typename intT, typename rationalT>
    struct cut_info
    {
        SINGLE_NUMERIC_TYPES_DEFS(intT, rationalT)

        is_vec ray;
        int_type rhs;
        rational_type frac_rhs;

        bool operator<(const cut_info &arg)
        {
            bool is_rays_less = Utility::SparseVecLess<int_type>()(ray, arg.ray);
            if (is_rays_less)
                return true;

            bool is_rays_equal = Utility::SparseVecEquiv<int_type>()(ray, arg.ray);
            if (is_rays_equal)
            {
                return rhs < arg.rhs;
            }

            return false;
        }

        bool operator==(const cut_info &arg)
        {
            return Utility::SparseVecEquiv<int_type>()(ray, arg.ray) && rhs == arg.rhs;
        }
    };
    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class cuts_gen_base : public dual_type_base_alg<intT, rationalT, lintT, lrationalT>
    {
    private:
        using my_type = cuts_gen_base<intT, rationalT, lintT, lrationalT>;
        using base = dual_type_base_alg<intT, rationalT, lintT, lrationalT>;

    public:
        DUAL_TYPE_ALGO_DEFS

        using cut_info = cut_info<intT, rationalT>;
        using local_instance = local_instance<intT, rationalT>;

    protected:
        std::vector<cut_info> _resulting_cuts;
        local_instance _instance_copy;

        vector<cut_info> generate_cuts_from_row_matrix(const is_mat &cuts_mat)
        {
            vector<cut_info> res;

            rs_vec rat_b = cuts_mat.template cast<rational_type>() * _instance_copy.v;
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
            out() << "Transforming rays to cuts..." << endl;
            timer.start();

            for (auto &mat : mat_pool)
            {
                mat = mat.transpose().eval();
            }

            for (const auto &mat : mat_pool)
            {
                auto new_cuts = generate_cuts_from_row_matrix(mat, _instance_copy.v);

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
        logger::stream &out()
        {
            log.out() << "ILPCuts-CutsGenerator: ";
            return log.out();
        }

        void set_instance(local_instance ins) { _instance_copy = std::move(ins); }

        virtual void generate_cuts() = 0;

        const std::vector<cut_info> &get_cuts() const { return _resulting_cuts; }

        // this function takes a collection of resulting cuts
        // and prints them to a string
        // each cut is separtated by \n
        // the last cut is also finished by \n
        void export_cuts(ofstream &out_file)
        {
            for (size_t i = 0; i < _resulting_cuts.size(); ++i)
            {
                out_file << Utility::sparse_vec_to_str(_resulting_cuts[i].ray) << " "
                         << _resulting_cuts[i].rhs << "\n";
            }
        }

        void export_cuts_full(ofstream &out_file)
        {
            for (size_t i = 0; i < _resulting_cuts.size(); ++i)
            {
                out_file << Utility::dense_vec_to_str(_resulting_cuts[i].ray) << " "
                         << _resulting_cuts[i].rhs << " "
                         << _resulting_cuts[i].frac_rhs << "\n";
            }
        }

        void import_cuts(istream &in_file)
        {
            clear_cuts();

            string line;

            while (std::getline(in_file, line))
            {
                istringstream iss_line(line);
                vector<int_type> elems;
                while (!iss_line.eof())
                {
                    elems.emplace_back(int_type(0));
                    iss_line >> elems.back();
                }

                size_t dim = elems.size() - 1;
                is_vec new_ray(dim);
                for (size_t i = 0; i < dim; ++i)
                {
                    if (elems[i] != 0)
                        new_ray.insertBack(i) = elems[i];
                }
                _resulting_cuts.emplace_back(cut_info{new_ray, elems.back(), 0.0});
            }
        }

        void clear_cuts()
        {
            _resulting_cuts.clear();
        }

        void check_cuts_validity()
        {
            size_t n = _AT_copy.cols();
            assert(n > 0);

            for (size_t i = 0; i < _resulting_cuts.size(); ++i)
            {
                const auto &cut = _resulting_cuts[i];

                rs_vec t = Utility::system_solve(_instance_copy.A.transpose().template cast<rational_type>().eval(), cut.ray.template cast<rational_type>().eval());

                for (typename rs_vec::InnerIterator it(t); it; ++it)
                {
                    if (it.value() < 0 && abs(it.value()) > tolerance)
                    {
                        stringstream msg;
                        msg << "ILPCuts error: cut(" << i << ") is invalid: negative solution(" << it.value() << ")...\n";
                        throw Exception(msg.str());
                    }
                }

                // rational_type vall = cut.ray.template cast<rational_type>().dot(_optLP_copy);
                // out() << "vall: " << vall << endl;

                if (abs(cut.ray.template cast<rational_type>().dot(_instance_copy.v) - cut.frac_rhs) > err_tolerance)
                {
                    stringstream msg;
                    msg << "ILPCuts error: cut(" << i << ") is invalid: invalid fractional rhs...\n";
                    throw Exception(msg.str());
                }

                if (floor(cut.frac_rhs) != cut.rhs)
                {
                    stringstream msg;
                    msg << "ILPCuts error: cut(" << i << ") is invalid: invalid floor rounding...\n";
                    throw Exception(msg.str());
                }
            }
        }
    };
} // namespace ILPCuts
