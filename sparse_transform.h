#pragma once

#include "basic_utils.h"

namespace ILPCuts
{
    namespace Utility
    {
        template <typename intT>
        struct sparse_transform
        {
            using int_type = intT;

            size_t i;
            SparseVector<int_type> t;
            int_type denum;
        };

        template <typename int_type, typename rational_type>
        SparseVector<rational_type> compute_transformed_coeffs(const SparseVector<rational_type> &x,
                                                               const sparse_transform<int_type> &trans, double tolerance)
        {

            using rs_vec = SparseVector<rational_type>;
            using is_vec = SparseVector<int_type>;

            size_t n = x.rows();

            /*	cout << "N:   " << n << endl;
                cin.ignore();*/

            rs_vec res(n);

            int_type ti = trans.t.coeff(trans.i);
            rational_type xi = x.coeff(trans.i);

            if (xi == 0)
            {
                rs_vec xx = x * ti;
                res = xx / trans.denum;
            }
            else
            {
                typename rs_vec::InnerIterator x_it(x);
                typename is_vec::InnerIterator t_it(trans.t);
                bool negative_solution_found = false;

                while (t_it && x_it && !negative_solution_found)
                {
                    if (t_it.index() < x_it.index())
                    {
                        negative_solution_found = true;
                    }
                    else if (x_it.index() < t_it.index())
                    {
                        rational_type y(x_it.value() * ti / trans.denum);
                        res.insert(x_it.index()) = std::move(y);
                        ++x_it;
                    }
                    else
                    {
                        rational_type y((x_it.value() * ti - static_cast<rational_type>(t_it.value()) * xi) / trans.denum);
                        if (y < rational_type(-tolerance))
                            negative_solution_found = true;
                        else
                            res.insert(t_it.index()) = std::move(y);
                        /*if (y > tolerance)
                            res.insert(t_it.index()) = std::move(y);
                        else if (y < -tolerance)
                            negative_solution_found = true;*/
                        ++x_it;
                        ++t_it;
                    }
                }

                if (!negative_solution_found)
                {
                    if (t_it)
                    {
                        negative_solution_found = true;
                    }
                    else
                    {
                        while (x_it)
                        {
                            rational_type y(x_it.value() * ti / trans.denum);
                            res.insert(x_it.index()) = std::move(y);
                            ++x_it;
                        }

                        res.coeffRef(trans.i) = xi;
                    }
                }

                if (negative_solution_found)
                    res.resize(0);
            }

            return std::move(res);
        }

        template <typename int_type>
        SparseMatrix<int_type> compute_adj_Q(size_t n, const sparse_transform<int_type> &trans)
        {
            using is_mat = SparseMatrix<int_type>;
            using is_vec = SparseVector<int_type>;

            const auto &t = trans.t;
            size_t i = trans.i;
            const auto &denum = trans.denum;
            const auto &ti = t.coeff(i);

            is_mat MM(n, n);
            using Tripl = Triplet<int_type>;
            std::vector<Tripl> triplets;
            for (size_t j = 0; j < n; ++j)
            {
                if (j != i)
                    triplets.push_back(Tripl(j, j, ti));
                else
                    triplets.push_back(Tripl(j, j, denum));
            }

            for (typename is_vec::InnerIterator it(t); it; ++it)
            {
                if (it.value() > 0 && it.index() != i)
                    triplets.emplace_back(Tripl(
                        it.index(),
                        i,
                        int_type(-it.value())));
            }

            MM.setFromTriplets(triplets.begin(), triplets.end());

            return std::move(MM);
        }
    } // namespace Utility

    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class transform_base
    {
    public:
        using rational_type = rationalT;
        using int_type = intT;

        using lrational_type = lrationalT;
        using lint_type = lintT;

        using is_mat = SparseMatrix<int_type>;
        using rs_mat = SparseMatrix<rational_type>;

        using lis_mat = SparseMatrix<lint_type>;
        using lrs_mat = SparseMatrix<lrational_type>;

        using is_vec = SparseVector<int_type>;
        using rs_vec = SparseVector<rational_type>;

        using lis_vec = SparseVector<lint_type>;
        using lrs_vec = SparseVector<lrational_type>;

    protected:
        Utility::RatToInt_Unary<lint_type, lrational_type> _lrat_to_lint;
        Utility::RatToInt_Unary<int_type, rational_type> _rat_to_int;

        const double &_tolerance;

    public:
        transform_base(const double &a_tol, const double &a_err_tol)
            : _tolerance(a_tol), _lrat_to_lint(a_err_tol), _rat_to_int(a_err_tol) {}

        const double &get_tolerance() const { return _tolerance; }

        virtual is_vec operator()(const is_mat &A) const = 0;
        virtual lint_type operator[](const lint_type &det) const = 0;
        virtual rs_vec compute_coeffs(const rs_vec &old_coefs) const = 0;
        virtual size_t get_i() const = 0;
        virtual lis_mat compute_adj(lis_mat old_adjA) const = 0;
    };

    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    class sparse_transform : public transform_base<intT, rationalT, lintT, lrationalT>
    {
    private:
        using base = transform_base<intT, rationalT, lintT, lrationalT>;

        using typename base::int_type;
        using typename base::lint_type;
        using typename base::lrational_type;
        using typename base::rational_type;

        using typename base::is_mat;
        using typename base::lis_mat;
        using typename base::lrs_mat;
        using typename base::rs_mat;

        using typename base::is_vec;
        using typename base::lis_vec;
        using typename base::lrs_vec;
        using typename base::rs_vec;

        using base::_lrat_to_lint;
        using base::_rat_to_int;
        using base::_tolerance;

        Utility::sparse_transform<int_type> data;

    public:
        sparse_transform(const double &a_tol, const double &a_err_tol, size_t ai,
                         const is_vec &at, const int_type &adenum)
            : base(a_tol, a_err_tol), data{ai, at, adenum}
        {
        }

        is_vec operator()(const is_mat &A) const override
        {
            rs_vec res = (A * data.t).template cast<rational_type>() / data.denum;

            return res.unaryExpr(_rat_to_int);
        }

        lint_type operator[](const lint_type &det) const override
        {
            const int_type &ti = data.t.coeff(data.i);
            lrational_type res(static_cast<lrational_type>(det * ti) / data.denum);

            return _lrat_to_lint(res);
        }

        rs_vec compute_coeffs(const rs_vec &old_coefs) const override
        {
            return Utility::compute_transformed_coeffs(old_coefs,
                                                       data,
                                                       _tolerance);
        }

        size_t get_i() const override { return data.i; }

        lis_mat compute_adj(lis_mat old_adjA) const override
        {
            size_t n = old_adjA.cols();
            is_mat Q = compute_adj_Q(n, data);
            old_adjA = ((Q.template cast<lint_type>() * old_adjA).template cast<lrational_type>() / data.denum)
                           .unaryExpr(_lrat_to_lint)
                           .eval();

            return std::move(old_adjA);
        }
    };

    template <typename intT, typename rationalT, typename lintT, typename lrationalT>
    struct lsparse_transform : public transform_base<intT, rationalT, lintT, lrationalT>
    {
    private:
        using base = transform_base<intT, rationalT, lintT, lrationalT>;

        using typename base::int_type;
        using typename base::lint_type;
        using typename base::lrational_type;
        using typename base::rational_type;

        using typename base::is_mat;
        using typename base::lis_mat;
        using typename base::lrs_mat;
        using typename base::rs_mat;

        using typename base::is_vec;
        using typename base::lis_vec;
        using typename base::lrs_vec;
        using typename base::rs_vec;

        using base::_lrat_to_lint;
        using base::_rat_to_int;
        using base::_tolerance;

        Utility::sparse_transform<lint_type> data;

        static constexpr auto _to_single_int = conv<lint_type, int_type>::make;
        static constexpr auto _to_single_rat = conv<lrational_type, rational_type>::make;

        // Utility::ToSingle<lint_type, int_type> _to_single_int;
        // Utility::ToSingle<lrational_type, rational_type> _to_single_rat;

    public:
        lsparse_transform(const double &a_tol, const double &a_err_tol, size_t ai,
                          const lis_vec &at, const lint_type &adenum)
            : base(a_tol, a_err_tol), data{ai, at, adenum}
        {
        }

        is_vec operator()(const is_mat &A) const override
        {
            lrs_vec res = (A.template cast<lint_type>() * data.t).template cast<lrational_type>() / data.denum;

            return res.unaryExpr(_lrat_to_lint).unaryExpr(_to_single_int);
        }

        lint_type operator[](const lint_type &det) const override
        {
            const lint_type &ti = data.t.coeff(data.i);
            lrational_type res(static_cast<lrational_type>(det * ti) / data.denum);
            return _lrat_to_lint(res);
        }

        rs_vec compute_coeffs(const rs_vec &old_coefs) const override
        {
            return Utility::compute_transformed_coeffs(old_coefs.template cast<lrational_type>().eval(),
                                                       data,
                                                       _tolerance)
                .unaryExpr(_to_single_rat);
        }

        size_t get_i() const override { return data.i; }

        lis_mat compute_adj(lis_mat old_adjA) const override
        {
            size_t n = old_adjA.cols();
            lis_mat Q = compute_adj_Q(n, data);
            old_adjA = ((Q * old_adjA).template cast<lrational_type>() / data.denum)
                           .unaryExpr(_lrat_to_lint)
                           .eval();

            return std::move(old_adjA);
        }
    };

} // namespace ILPCuts
