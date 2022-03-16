#pragma once

#include"conic_cuts_base.h"

#include <omp.h>
#include <atomic>

namespace ILPCuts {

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	class ConicCutsGenerator :public ConicCutsGenerator_Base<intT,rationalT,lintT,lratT> {
	private:
		using base = typename ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>;

	public:
		using typename base::rational_type;
		using typename base::int_type;
		using typename base::lrational_type;
		using typename base::lint_type;

		using typename base::is_mat;
		using typename base::rs_mat;
		using typename base::lis_mat;
		using typename base::lrs_mat;

		using typename base::is_vec;
		using typename base::rs_vec;
		using typename base::lis_vec;
		using typename base::lrs_vec;

		using typename base::sparce_transform;
		using typename base::dir_coeffs;
		using typename base::cone_traits;
		using typename base::process_func;

		using base::tolerance;
		using base::threads;
		using base::out;
		using base::err;
		using base::err_tolerance;
		using base::cone_limit_per_direction;
		using base::direction_strategies;
		using base::_rat_to_int;
		using base::_lrat_to_lint;
		using base::rgen;
		using base::_to_single_int;
		using base::_to_single_rat;

		using base::base;

	private:

		struct cone_traits_spec
			:public base::cone_traits {
			lis_mat adjA;
		};

		std::unique_ptr<cone_traits> create_cone_traits(const is_mat& A, const lis_mat& adjA, const lint_type& det,
			const std::vector<rs_vec>& directions_pool) {
			auto res = std::make_unique<cone_traits_spec>();

			res -> A = A;
			res -> det = det;
			res -> adjA = adjA;

			for (size_t i = 0; i < directions_pool.size(); ++i) {
				res -> dir_coeffs_pull.emplace_back(dir_coeffs{
					(adjA.cast<lrational_type>() * directions_pool[i]).unaryExpr(_to_single_rat), i
					});
			}

			return std::move(res);
		}

		//we can run this method only if dir_coeffs_pull is not empty
		std::unique_ptr<cone_traits> compute_cone_traits(const cone_traits * base_C,
			const sparce_transform& trans,
			std::vector<dir_coeffs> dir_coeffs_pull) override {

			//++directions_counter[dir_coeffs_pull[0].dir_ID];

			const cone_traits_spec* C = static_cast<const cone_traits_spec*>(base_C);

			size_t n = C->A.cols();
			const auto& t = trans.t;
			size_t i = trans.i;
			const auto& denum = trans.denum;

			auto new_cone = std::make_unique<cone_traits_spec>();

			new_cone -> dir_coeffs_pull = std::move(dir_coeffs_pull);

			try {
				new_cone -> A = C -> A;
				/*new_cone.A.col(i) = ((C.A * t).cast<rational_type>() / denum)
					.unaryExpr(*_rat_to_int);*/
				new_cone -> A.col(i) = ( (C->A.cast<lint_type>() * t).cast<lrational_type>() / denum)
					.unaryExpr(*_lrat_to_lint).unaryExpr(_to_single_int);
				Utility::prune_zeroes_int(new_cone -> A);
			}
			catch (ToleranceErrorException exc) {
				out() << "A-column construction tolerance error: " << exc.what() << '\n';
				err() << "A-column construction tolerance error: " << exc.what() << '\n';
			}

			const auto & ti = t.coeff(i);
			try {
				new_cone -> det = (*_lrat_to_lint)(static_cast<lrational_type>(C->det * ti) / denum);
			}
			catch (ToleranceErrorException exc) {
				out() << "Determinant recompute tolerance error: " << exc.what() << '\n';
				err() << "Determinant recompute tolerance error: " << exc.what() << '\n';
			}

			/*if (new_cone.det == 0) {
				out() << det <<", "<< ti << std::endl;
				std::cin.ignore();
			}*/

			lis_mat MM(n, n);
			using Tripl = Triplet<lint_type>;
			std::vector<Tripl> triplets;
			for (size_t j = 0; j < n; ++j) {
				if (j != i)
					triplets.push_back(Tripl(j, j, ti));
				else
					triplets.push_back(Tripl(j, j, denum));
			}

			for (typename lis_vec::InnerIterator it(t); it; ++it) {
				if (it.value() > 0 && it.index() != i)
					triplets.push_back(Tripl(it.index(), i, -it.value()));
			}

			MM.setFromTriplets(triplets.begin(), triplets.end());

			try {
				new_cone -> adjA = ((MM * C->adjA).cast<lrational_type>() / denum).unaryExpr(*_lrat_to_lint);
				Utility::prune_zeroes_int(new_cone -> adjA);
			}
			catch (ToleranceErrorException exc) {
				out() << "Adj. recompute tolerance error: " << exc.what() << '\n';
				err() << "Adj. recompute tolerance error: " << exc.what() << '\n';
			}

			return std::move(new_cone);
		}

		void process_node(const cone_traits& base_C, process_func pf) override {

			const cone_traits_spec& C = static_cast<const cone_traits_spec&>(base_C);

			const is_mat& A = C.A;
			const lint_type& det = C.det;
			const lis_mat& adjA = C.adjA;
			size_t n = A.cols();

			if (det % 2 == 0) {
				lis_vec pivot;
				for (size_t i = 0; i < n; ++i) {
					pivot = adjA.col(i);
					lint_type adj_col_2pow = Utility::vec_max2pow(pivot);
					if ((static_cast<lint_type>(det) / adj_col_2pow) % 2 == 0) {
						pivot /= adj_col_2pow;
						break;
					}
				}

				pivot = pivot.unaryExpr(std::ref(Utility::mod2_unaryop<lint_type>)).eval();

				/*for (size_t i = 0; i < pivot.size(); ++i)
					if (pivot[i] != 0 && pivot[i] != 1) {
						cout << "pivot[i]" << pivot[i] << endl;
						cin.ignore();
					}*/

				for (typename is_vec::InnerIterator it(pivot); it; ++it) {
					if (it.value() > 0)
						pf(sparce_transform{
							static_cast<size_t>(it.index()), std::move(pivot), 2
							});
				}

			}
			else {
				//std::cout << "oddddddddddddddddodddd" << std::endl;

				lis_vec pivot;
				//out() << A * adjA << std::endl;
				for (size_t i = 0; i < n; ++i) {
					pivot = adjA.col(i).unaryExpr(Utility::mod_unaryop(det));
					bool is_any_nonzero = false;
					for (typename lis_vec::InnerIterator it(pivot); it; ++it) {
						if (it.value() > 0) {
							is_any_nonzero = true;
							break;
						}
					}
					//out() << i << " ";

					if (is_any_nonzero)
						break;
				}

				for (typename lis_vec::InnerIterator it(pivot); it; ++it) {
					if (it.value() > 0) {
						lint_type lambda = Utility::mod_equation_solve<lint_type>(it.value(), det - gcd(det, it.value()), det);
						lis_vec t = (pivot * lambda).unaryExpr(Utility::mod_unaryop(det));

						pf(sparce_transform{
							static_cast<size_t>(it.index()), std::move(t), std::move(det)
							});
					}
				}
			}
		}

	public:

		/*ConicCutsGenerator(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:base(outs, errs) {}*/

		std::vector<is_vec> generate_cuts(const is_mat& A, const lis_mat& adjA, const rs_vec& optLP, const lint_type& det);

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP) override;
	};

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT, rationalT, lintT, lratT>::is_vec>
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat& A, const lis_mat& adjA, const rs_vec& optLP, const lint_type& det) {

		auto directions_pool = std::move(base::generate_directions_pool(A));

		return std::move(base::generate_cuts( 
			create_cone_traits(A, adjA, det, directions_pool),
			optLP
		));
	}

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT, rationalT, lintT, lratT>::is_vec>
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat& A, const rs_vec& optLP) {

		lrational_type det = 0;
		lrs_mat adjA;

		try {
			auto start = std::chrono::steady_clock::now();
			std::tie(adjA, det) = std::move(
				Utility::compute_abs_adjugate_and_det(A.cast<lrational_type>().eval())
			);
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			out() << "ConicCuts: matrix inversion time: " << elapsed_seconds.count() << "s\n";

			return std::move(generate_cuts(A, adjA.unaryExpr(*_lrat_to_lint), optLP, (*_lrat_to_lint)(det) ));
		}
		catch (ToleranceErrorException exc) {
			out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
		}
	}
}