#pragma once

#include "base_algo.h"
#include "local_instance.h"
#include "Eigen_solving_utilities.h"

namespace ILPCuts
{
    struct cut_info
    {
        rs_vec ray;
        rational_type rhs;
        rational_type precut_rhs;
        std::vector<size_t> intComp_idxs;

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
    class cuts_gen_base : public base_alg
    {
    public:
        local_instance LocalInstance;

    protected:
        std::vector<cut_info> _resulting_cuts;

    public:
        cuts_gen_base(logger &a_log) : base_alg(a_log) {}

        logger::stream &out()
        {
            log.out() << "ILPCuts-CutsGenerator: ";
            return log.out();
        }

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
                out_file << Utility::sparse_vec_to_str(_resulting_cuts[i].ray) << " "
                         << _resulting_cuts[i].rhs << " "
                         << _resulting_cuts[i].precut_rhs << "\n"
                         << Utility::dense_vec_to_str(_resulting_cuts[i].intComp_idxs) << "\n\n";
            }
        }

        void import_cuts(istream &in_file)
        {
            clear_cuts();

            string line;

            while (std::getline(in_file, line))
            {
                istringstream iss_line(line);
                vector<rational_type> elems;
                while (!iss_line.eof())
                {
                    elems.emplace_back(rational_type(0.0));
                    iss_line >> elems.back();
                }

                size_t dim = elems.size() - 1;
                rs_vec new_ray(dim);
                for (size_t i = 0; i < dim; ++i)
                {
                    if (is_zero(elems[i], tolerance))
                        new_ray.insertBack(i) = elems[i];
                }
                _resulting_cuts.emplace_back(cut_info{new_ray, elems.back(), 0.0, std::vector<size_t>()});
            }
        }

        // rs_mat adjust_normal_cone(rs_mat A, const rs_vec &c)
        // {
        //     rs_vec t = Utility::system_solve(A, c);

        //     for (typename rs_vec::InnerIterator it(t); it; ++it)
        //     {
        //         if (it.value() < 0 && !is_zero(it.value(), tolerance))
        //         {
        //             A.col(it.index()) = -A.col(it.index());
        //         }
        //     }

        //     return A;
        // }

        void clear_cuts()
        {
            _resulting_cuts.clear();
        }

        void check_pure_int_cuts_validity()
        {
            is_mat AT = LocalInstance.A.transpose().unaryExpr(rat_to_int);

            size_t n = AT.cols();
            assert(n > 0);

            for (size_t i = 0; i < _resulting_cuts.size(); ++i)
            {
                const auto &cut = _resulting_cuts[i];

                rs_vec t = Utility::system_solve(AT.cast<rational_type>().eval(), cut.ray.cast<rational_type>().eval());

                for (typename rs_vec::InnerIterator it(t); it; ++it)
                {
                    if (it.value() < 0 && abs(it.value()) > tolerance)
                    {
                        stringstream msg;
                        msg << "ILPCuts error: cut(" << i << ") is invalid: negative solution(" << it.value() << ")...\n";
                        throw Exception(msg.str());
                    }
                }

                if (abs(cut.ray.cast<rational_type>().dot(LocalInstance.v) - cut.precut_rhs) > err_tolerance)
                {
                    stringstream msg;
                    msg << "ILPCuts error: cut(" << i << ") is invalid: invalid fractional rhs...\n";
                    throw Exception(msg.str());
                }

                if (floor(cut.precut_rhs) != cut.rhs)
                {
                    stringstream msg;
                    msg << "ILPCuts error: cut(" << i << ") is invalid: invalid floor rounding...\n";
                    throw Exception(msg.str());
                }
            }
        }
    };
} // namespace ILPCuts
