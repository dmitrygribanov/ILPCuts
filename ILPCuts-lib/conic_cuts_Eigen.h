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

		/*using id_mat = Matrix<int_type, Dynamic, Dynamic>;
		using rd_mat = Matrix<rational_type, Dynamic, Dynamic>;

		using id_vec = Matrix<int_type, Dynamic, 1>;
		using rd_vec = Matrix<rational_type, Dynamic, 1>;*/

	private:

		std::vector<is_mat> partial_uni_decomp
		(const is_mat& A, const std::vector<rs_vec>& directions_pool,
			const lis_mat& adjA, const lint_type& det);

		struct cone_traits
			:public base::cone_traits {
			lis_mat adjA;
		};

		cone_traits create_cone_traits(const is_mat& A, const lis_mat& adjA, const lint_type& det,
			const std::vector<rs_vec>& directions_pool) {
			cone_traits res;

			res.A = A;
			res.det = det;
			res.adjA = adjA;

			for (size_t i = 0; i < directions_pool.size(); ++i) {
				res.dir_coeffs_pull.emplace_back(dir_coeffs{
					adjA.cast<rational_type>() * directions_pool[i], i
					});
			}

			return std::move(res);
		}

		/*cone_traits create_cone_traits(const is_mat& A, const std::vector<rs_vec>& directions_pool) {
			is_mat adjA;
			int_type det;
			std::tie(adjA, det) = std::move(
				Utility::compute_abs_adjugate_and_det(A.cast<rational_type>().eval())
			);
			Utility::prune_zeroes_int(adjA);

			return std::move(create_cone_traits(A, adjA, det, directions_pool));
		}*/

		//we can run this method only if dir_coeffs_pull is not empty
		cone_traits compute_cone_traits(const cone_traits& C,
			const sparce_transform& trans,
			//std::vector<size_t> & directions_counter,
			std::vector<dir_coeffs> dir_coeffs_pull) {

			//++directions_counter[dir_coeffs_pull[0].dir_ID];

			size_t n = C.A.cols();
			const auto& t = trans.t;
			size_t i = trans.i;
			const auto& denum = trans.denum;

			cone_traits new_cone;

			new_cone.dir_coeffs_pull = std::move(dir_coeffs_pull);

			try {
				new_cone.A = C.A;
				/*new_cone.A.col(i) = ((C.A * t).cast<rational_type>() / denum)
					.unaryExpr(*_rat_to_int);*/
				new_cone.A.col(i) = ( (C.A.cast<lint_type>() * t).cast<lrational_type>() / denum)
					.unaryExpr(*_lrat_to_lint).unaryExpr(_to_single_int);
				Utility::prune_zeroes_int(new_cone.A);
			}
			catch (ToleranceErrorException exc) {
				out() << "A-column construction tolerance error: " << exc.what() << '\n';
				err() << "A-column construction tolerance error: " << exc.what() << '\n';
			}

			int_type ti = t.coeff(i);
			try {
				new_cone.det = (*_lrat_to_lint)(static_cast<lrational_type>(C.det * ti) / denum);
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
				new_cone.adjA = ((MM * C.adjA).cast<lrational_type>() / denum).unaryExpr(*_lrat_to_lint);
				Utility::prune_zeroes_int(new_cone.adjA);
			}
			catch (ToleranceErrorException exc) {
				out() << "Adj. recompute tolerance error: " << exc.what() << '\n';
				err() << "Adj. recompute tolerance error: " << exc.what() << '\n';
			}

			return std::move(new_cone);
		}

		std::vector<cone_traits> compute_layer(const cone_traits& cone, const std::vector<rs_vec>& directions_pool, 
			std::vector<Atomic<size_t>>& directions_evaluation_counter);

	public:

		/*ConicCutsGenerator(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:base(outs, errs) {}*/

		std::vector<is_vec> generate_cuts(const is_mat& A, const lis_mat& adjA, const rs_vec& optLP, const lint_type& det);

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP) override;
	};

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT,rationalT,lintT,lratT>::is_mat> 
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::partial_uni_decomp
		(const is_mat& A, const std::vector<rs_vec>& directions_pool,
			const lis_mat& adjA, const lint_type& det) {
		//Задаю глобальную параллельность openmp
		omp_set_dynamic(max(threads, size_t(1)));
		std::vector<cone_traits> current_layer;
		std::vector<std::vector<is_mat>> result(threads, std::vector<is_mat>());

		current_layer.push_back(std::move(create_cone_traits(
			A, adjA, det, directions_pool
		) ));

		while (!current_layer.empty()) {

			/*out() << "Curr level size: " << current_layer.size() << std::endl;
			std::cin.ignore();*/

			std::vector<Atomic<size_t>> directions_evaluation_counter(directions_pool.size(), Atomic<size_t>(0));

			int dynamic_threads = std::min(threads, current_layer.size());

			std::vector<std::vector<cone_traits>> next_layers(dynamic_threads, std::vector<cone_traits>());
			omp_set_num_threads(dynamic_threads);
#ifndef ILPCUTS_PARALLEL_OFF
#pragma omp parallel for
#endif
			for (int i = 0; i < current_layer.size(); ++i) {
				const cone_traits& node = current_layer[i];
				std::vector<cone_traits> next_node_layer = std::move(compute_layer(node, directions_pool, directions_evaluation_counter));
				int tid = omp_get_thread_num();
				for (cone_traits& next_node : next_node_layer) {
					if (next_node.det == 1)
						result[tid].push_back(std::move(next_node.A));
					else
						next_layers[tid].push_back(std::move(next_node));
				}
			}

			std::swap(current_layer, next_layers[0]);
			for (int i = 1; i < next_layers.size(); ++i) {
				current_layer.insert(
					current_layer.end(),
					std::make_move_iterator(next_layers[i].begin()),
					std::make_move_iterator(next_layers[i].end())
				);
			}
		}

		std::vector<is_mat>& concatenated_result = result[0];
		for (int i = 1; i < result.size(); ++i) {
			concatenated_result.insert(
				concatenated_result.end(),
				std::make_move_iterator(result[i].begin()),
				std::make_move_iterator(result[i].end())
			);
		}

		return std::move(concatenated_result);
	}

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT, rationalT, lintT, lratT>::is_vec>
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat& A, const lis_mat& adjA, const rs_vec& optLP, const lint_type& det) {

		auto start = std::chrono::steady_clock::now();

		std::vector<is_vec> res_cuts;

		std::vector<rs_vec> directions_pool;
		for (auto & strategy : direction_strategies) {
			std::vector<rs_vec> strategy_directions = std::move(strategy->generate_directons(A));
			//out() << strategy_directions[0] << std::endl;
			directions_pool.insert(directions_pool.end(), strategy_directions.begin(), strategy_directions.end());
		}

		//auto directions_pool = std::move(directions_strategy.generate_directons(A));

		std::vector<is_mat> uni_mat_pool = std::move(partial_uni_decomp(A, directions_pool,
			adjA, det));

		for (is_mat& uni_mat : uni_mat_pool) {
			uni_mat = uni_mat.transpose();
		}
		

		for (const is_mat& uni_mat : uni_mat_pool) {
			rs_vec rat_b = uni_mat.cast<rational_type>() * optLP;
			is_vec rounded_b = rat_b.unaryExpr(Utility::Round<rational_type>()).cast<int_type>();
			is_vec floor_b = rat_b.unaryExpr(Utility::Floor<rational_type>()).cast<int_type>();

			rs_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).cwiseAbs();

			/*std::cout << uni_mat << std::endl;
			std::cout << optLP.transpose() << std::endl;
			std::cout << rat_b.transpose() << std::endl;
			std::cout << rounded_b.transpose() << std::endl;
			std::cout << floor_b.transpose() << std::endl;
			std::cout << diff_b.transpose() << std::endl;*/

			typename rs_vec::InnerIterator diff_it(diff_b);
			//rs_vec::InnerIterator floor_it(floor_b);
			//size_t jjjj = 0;
			while(diff_it) {
				/*if (diff_it.value() != 0)
				{
					cout << ++jjjj << ", " << rat_b.coeff(diff_it.index()) - floor_b.coeff(diff_it.index()) << endl;
					cin.ignore();
				}*/
				if (diff_it.value() > tolerance) {
					is_vec cut = uni_mat.row(diff_it.index());
					//cut.reserve(1);
					cut.insert(cut.size()) = floor_b.coeff(diff_it.index());
					//cout << cut << endl;
					//out() << "cut has been added: (" << cut.transpose() << ")\n";
					res_cuts.push_back(std::move(cut));
				}
				++diff_it;
				//++floor_it;
			}

			/*for (size_t i = 0; i < diff_b.size(); ++i) {
				if (diff_b[i] > tolerance) {
					is_vec cut = uni_mat.row(i);
					cut.conservativeResize(cut.size() + 1);
					cut[cut.size() - 1] = floor_b[i];
					out() << "cut has been added: (" << cut.transpose() << ")\n";
					res_cuts.push_back(std::move(cut));
				}
			}*/
		}

		std::sort(res_cuts.begin(), res_cuts.end(), Utility::SparseVecComp<rs_vec>());
		auto last = std::unique(res_cuts.begin(), res_cuts.end(), Utility::SparseVecEquiv<rs_vec>());
		// v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
		res_cuts.erase(last, res_cuts.end());

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		out() << "ConicCuts: cuts generation time: " << elapsed_seconds.count() << "s\n";
		out() << "ConicCuts: total cuts: " << res_cuts.size() << "\n";

		return std::move(res_cuts);
	}

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT, rationalT, lintT, lratT>::is_vec>
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::generate_cuts(const is_mat& A, const rs_vec& optLP) {

		rational_type det = 0;
		rs_mat adjA;

		try {
			auto start = std::chrono::steady_clock::now();
			std::tie(adjA, det) = std::move(
				Utility::compute_abs_adjugate_and_det<rational_type>(A.cast<rational_type>())
			);
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			out() << "ConicCuts: matrix inversion time: " << elapsed_seconds.count() << "s\n";

			return std::move(generate_cuts(A, adjA.unaryExpr(*_rat_to_int), optLP, (*_rat_to_int)(det) ));
		}
		catch (ToleranceErrorException exc) {
			out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
		}
	}

	template<typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator<intT, rationalT, lintT, lratT>::cone_traits>
		ConicCutsGenerator<intT, rationalT, lintT, lratT>::compute_layer(const cone_traits& cone, const std::vector<rs_vec>& directions_pool,
			std::vector<Atomic<size_t>> & directions_evaluation_counter) {

		std::vector<cone_traits> result;

		const is_mat& A = cone.A;
		const is_mat& adjA = cone.adjA;
		const int_type& det = cone.det;
		size_t n = A.cols();

		out() << "Current determinant: " << det << std::endl;
		out() << "Nonzeros: " << A.nonZeros() << std::endl;

		auto process_new_candidate = [&](const sparce_transform& trans)
		{
			auto dirr_coeffs_pull = std::move(base::compute_dirr_coeffs_pull(
				cone, trans, directions_evaluation_counter
			));

			if (!dirr_coeffs_pull.empty()) {
				++directions_evaluation_counter[dirr_coeffs_pull[0].dir_ID];

				cone_traits new_cone = std::move(compute_cone_traits(
					cone, trans, dirr_coeffs_pull
				));

				result.push_back(std::move(new_cone));
				out() << "Adding step done..." << std::endl;
			}
		};

		if (det % 2 == 0) {
			//std::cout << "evennnnnnnnnnnnnnn" << std::endl;

			is_vec pivot;
			for (size_t i = 0; i < n; ++i) {
				pivot = adjA.col(i);
				int_type adj_col_2pow = Utility::vec_max2pow(pivot);
				if ((static_cast<int_type>(det) / adj_col_2pow) % 2 == 0) {
					pivot /= adj_col_2pow;
					break;
				}
			}

			pivot = pivot.unaryExpr(std::ref(Utility::mod2_unaryop<int_type>)).eval();

			/*for (size_t i = 0; i < pivot.size(); ++i)
				if (pivot[i] != 0 && pivot[i] != 1) {
					cout << "pivot[i]" << pivot[i] << endl;
					cin.ignore();
				}*/

			for (typename is_vec::InnerIterator it(pivot); it; ++it) {
				if (it.value() > 0)
					process_new_candidate(sparce_transform{ 
						static_cast<size_t>(it.index()), std::move(pivot), 2 
					});
			}

		}
		else {
			//std::cout << "oddddddddddddddddodddd" << std::endl;

			is_vec pivot;
			//out() << A * adjA << std::endl;
			for (size_t i = 0; i < n; ++i) {
				pivot = adjA.col(i).unaryExpr(Utility::mod_unaryop(det));
				bool is_any_nonzero = false;
				for (typename is_vec::InnerIterator it(pivot); it; ++it) {
					if (it.value() > 0) {
						is_any_nonzero = true;
						break;
					}
				}
				//out() << i << " ";
				
				if (is_any_nonzero)
					break;
			}

			for (typename is_vec::InnerIterator it(pivot); it; ++it) {
				if (it.value() > 0) {
					int_type lambda = Utility::mod_equation_solve<int_type>(it.value(), det - gcd(det, it.value()), det);
					is_vec t = (pivot * lambda).unaryExpr(Utility::mod_unaryop(det));

					process_new_candidate(sparce_transform{ 
						static_cast<size_t>(it.index()), std::move(t), std::move(det) 
					});
				}
			}
		}

		return std::move(result);
	}
}