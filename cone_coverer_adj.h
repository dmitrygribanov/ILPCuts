#pragma once

#include "cone_coverer_base.h"
#include "Eigen_solving_utilities.h"
#include "ntl_module.h"
namespace ILPCuts
{
	class ConeCoverer_adj : public ConeCoverer_Base
	{
	private:
		using base = ConeCoverer_Base;

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
					(adjA.template cast<lrational_type>() * directions_pool[i].template cast<lrational_type>()).unaryExpr(to_single_rat), i});
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
				out() << "Adj. matrix recompute tolerance error: " << exc.what() << endl;
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

			is_vec res = pivot.unaryExpr(std::ref(Utility::mod2_unaryop<lint_type>)).unaryExpr(to_single_int);
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
				out() << "even step...\n";
				process_even_case(C, pf);
			}
			else
			{
				out() << "odd step...\n";
				process_odd_case(C, pf);
			}
		}

	public:
		using base::base;
		/*ConeCoverer_adj(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:base(outs, errs) {}*/

		std::vector<is_mat> generate_cuts(const is_mat &A, const lis_mat &adjA, const lint_type &det);

		std::vector<is_mat> generate_cuts(const is_mat &A) override;
	};

	std::vector<is_mat> ConeCoverer_adj::generate_cuts(const is_mat &A, const lis_mat &adjA, const lint_type &det)
	{

		out() << "Generating directions pool ..." << endl;

		auto directions_pool = std::move(base::generate_directions_pool(A));

		return base::generate_cuts(create_cone_traits(A, adjA, det, directions_pool));
	}

	std::vector<is_mat> ConeCoverer_adj::generate_cuts(const is_mat &A)
	{
		try
		{
			Timer timer;
			lint_type det(0);
			lis_mat adjA;

			out() << "matrix inversion is started..." << endl;
			timer.start();
			std::tie(adjA, det) = ModuleNTL::compute_inverse(A.template cast<lint_type>().eval());
			// std::tie(adjA, det) = std::move(
			// 	Utility::compute_abs_adjugate_and_det(A.template cast<lrational_type>().eval()));
			out() << "matrix inversion time: " << timer.finish_str();

			return generate_cuts(A, adjA, abs(det));
		}
		catch (ToleranceErrorException exc)
		{
			out() << "determinant computing tolerance error: " << exc.what() << endl;
			return std::vector<is_mat>();
		}
	}
} // namespace ILPCuts