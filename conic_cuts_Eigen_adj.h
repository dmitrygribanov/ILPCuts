#pragma once

#include "conic_cuts_base.h"
#include "Eigen_solving_utilities.h"

#include <omp.h>
#include <atomic>

namespace ILPCuts
{

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	class ConicCutsGenerator_adj : public ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>
	{
	private:
		using base = ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>;

	public:
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

		using typename base::cone_traits;
		using typename base::dir_coeffs;
		using typename base::lsparse_transform;
		using typename base::process_func;
		using typename base::sparse_transform;
		using typename base::transform_base;

		using base::_lrat_to_lint;
		using base::_rat_to_int;
		using base::_to_single_int;
		using base::_to_single_rat;
		using base::cone_limit_per_direction;
		using base::direction_strategies;
		using base::err;
		using base::err_tolerance;
		using base::out;
		using base::rgen;
		using base::threads;
		using base::tolerance;

		using base::base;

	private:
		struct cone_traits_spec
			: public base::cone_traits
		{
			lis_mat adjA;
		};

		std::unique_ptr<cone_traits> create_cone_traits(const is_mat &A, const lis_mat &adjA, const lint_type &det,
														const std::vector<rs_vec> &directions_pool)
		{

			auto res = std::make_unique<cone_traits_spec>();

			res->A = A;
			res->det = det;
			res->adjA = adjA;

			for (size_t i = 0; i < directions_pool.size(); ++i)
			{
				res->dir_coeffs_pull.emplace_back(dir_coeffs{
					(adjA.template cast<lrational_type>() * directions_pool[i].template cast<lrational_type>()).unaryExpr(_to_single_rat), i});
			}

			return std::move(res);
		}

		// we can run this method only if dir_coeffs_pull is not empty
		std::unique_ptr<cone_traits> compute_cone_traits(const cone_traits *base_C,
														 const transform_base &trans,
														 std::vector<dir_coeffs> dir_coeffs_pull) override
		{

			const cone_traits_spec *C = static_cast<const cone_traits_spec *>(base_C);

			auto new_cone = std::make_unique<cone_traits_spec>();
			new_cone->move_from_base_ptr(base::compute_cone_traits(base_C, trans, dir_coeffs_pull));

			try
			{
				new_cone->adjA = trans.compute_adj(C->adjA);
				Utility::prune_zeroes_int(new_cone->adjA);
			}
			catch (ToleranceErrorException exc)
			{
				out() << "Adj. recompute tolerance error: " << exc.what() << '\n';
				err() << "Adj. recompute tolerance error: " << exc.what() << '\n';
			}

			return std::move(new_cone);
		}

		void process_even_case(const cone_traits_spec &C, process_func pf)
		{

			const is_mat &A = C.A;
			const lint_type &det = C.det;
			const lis_mat &adjA = C.adjA;
			size_t n = A.cols();

			lis_vec pivot;
			for (size_t i = 0; i < n; ++i)
			{
				pivot = adjA.col(i);
				lint_type adj_col_2pow = Utility::vec_max2pow(pivot);
				if ((static_cast<lint_type>(det) / adj_col_2pow) % lint_type(2) == 0)
				{
					pivot /= adj_col_2pow;
					break;
				}
			}

			is_vec res = pivot.unaryExpr(std::ref(Utility::mod2_unaryop<lint_type>)).unaryExpr(_to_single_int);
			Utility::prune_zeroes_int(res);

			for (typename is_vec::InnerIterator it(res); it; ++it)
			{
				if (it.value() > 0)
					pf(sparse_transform(tolerance, err_tolerance, static_cast<size_t>(it.index()), std::move(res), int_type(2)));
			}
		}

		void process_odd_case(const cone_traits_spec &C, process_func pf)
		{

			const is_mat &A = C.A;
			const lint_type &det = C.det;
			const lis_mat &adjA = C.adjA;
			size_t n = A.cols();

			lis_vec pivot;
			for (size_t i = 0; i < n; ++i)
			{
				pivot = adjA.col(i).unaryExpr(Utility::mod_unaryop(det));

				Utility::prune_zeroes_int(pivot);

				if (pivot.nonZeros() > 0)
					break;
			}

			for (typename lis_vec::InnerIterator it(pivot); it; ++it)
			{
				if (it.value() > 0)
				{
					lint_type lambda = Utility::mod_equation_solve<lint_type>(it.value(), lint_type(det - gcd(det, it.value())), det);
					lis_vec t = (pivot * lambda).unaryExpr(Utility::mod_unaryop(det));

					pf(lsparse_transform(tolerance, err_tolerance, static_cast<size_t>(it.index()), std::move(t), std::move(det)));
				}
			}
		}

		void process_node(const cone_traits &base_C, process_func pf) override
		{

			const cone_traits_spec &C = static_cast<const cone_traits_spec &>(base_C);

			if (C.det % lint_type(2) == 0)
			{
				out() << "even step..." << std::endl;
				process_even_case(C, pf);
			}
			else
			{
				out() << "odd step..." << std::endl;
				process_odd_case(C, pf);
			}
		}

	public:
		/*ConicCutsGenerator_adj(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:base(outs, errs) {}*/

		std::vector<is_vec> generate_cuts(const is_mat &A, const lis_mat &adjA, const rs_vec &optLP, const lint_type &det);

		std::vector<is_vec> generate_cuts(const is_mat &A, const rs_vec &optLP) override;
	};

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator_adj<intT, rationalT, lintT, lratT>::is_vec>
	ConicCutsGenerator_adj<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat &A, const lis_mat &adjA, const rs_vec &optLP, const lint_type &det)
	{

		auto directions_pool = std::move(base::generate_directions_pool(A));

		return std::move(base::generate_cuts(
			create_cone_traits(A, adjA, det, directions_pool),
			optLP));
	}

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator_adj<intT, rationalT, lintT, lratT>::is_vec>
	ConicCutsGenerator_adj<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat &A, const rs_vec &optLP)
	{
		try
		{
			lrational_type det(0.0);
			lrs_mat adjA;

			auto start = std::chrono::steady_clock::now();
			std::tie(adjA, det) = std::move(
				Utility::compute_abs_adjugate_and_det(A.template cast<lrational_type>().eval()));
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			out() << "ConicCuts: matrix inversion time: " << elapsed_seconds.count() << "s\n";

			return std::move(generate_cuts(A, adjA.unaryExpr(_lrat_to_lint), optLP, _lrat_to_lint(det)));
		}
		catch (ToleranceErrorException exc)
		{
			out() << "ILPCuts generator: determinant computing tolerance error: " << exc.what() << '\n';
			err() << "ILPCuts generator: determinant computing tolerance error: " << exc.what() << '\n';

			return std::vector<is_vec>();
		}
	}
}