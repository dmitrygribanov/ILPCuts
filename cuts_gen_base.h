#pragma once

#include "base_algo.h"
#include "local_instance.h"
#include "Eigen_solving_utilities.h"

namespace ILPCuts
{
    struct cut_info
    {
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

        // void set_instance(local_instance ins) { LocalInstance = std::move(ins); }

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

        is_mat adjust_normal_cone(is_mat A, const rs_vec &c)
        {
            rs_vec t = Utility::system_solve(A.template cast<rational_type>().eval(), c);

            for (typename rs_vec::InnerIterator it(t); it; ++it)
            {
                // it.value() < tolerance previouse code. working fast but its must be not correct
                if (it.value() < 0 && abs(it.value()) > tolerance)
                {
                    A.col(it.index()) = -A.col(it.index());
                }
            }

            return std::move(A);
        }

        void clear_cuts()
        {
            _resulting_cuts.clear();
        }

        void check_cuts_validity()
        {
            is_mat AT = LocalInstance.A.transpose();

            size_t n = AT.cols();
            assert(n > 0);

            for (size_t i = 0; i < _resulting_cuts.size(); ++i)
            {
                const auto &cut = _resulting_cuts[i];

                rs_vec t = Utility::system_solve(AT.template cast<rational_type>().eval(), cut.ray.template cast<rational_type>().eval());

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

                if (abs(cut.ray.template cast<rational_type>().dot(LocalInstance.v) - cut.frac_rhs) > err_tolerance)
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
