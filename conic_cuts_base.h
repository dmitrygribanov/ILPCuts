#pragma once

#include "base_module.h"
#include "Hilbert_base_alg.h"

// all optional modules are added at the moment
// to simplify experiments
#include "Eigen_solving_utilities.h"
#include "flint_module.h"
#include "ntl_module.h"

#include <omp.h>
#include <atomic>

namespace ILPCuts
{
	template <typename intT, typename rationalT, typename lintT, typename lrationalT>
	class ConicCutsGenerator_Base
		: public base_alg<intT, rationalT, lintT, lrationalT>
	{
	private:
		using my_type = ConicCutsGenerator_Base<intT, rationalT, lintT, lrationalT>;
		using base = base_alg<intT, rationalT, lintT, lrationalT>;

		friend class CPlexAdapter<my_type>;

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

		using base::_lrat_to_lint;
		using base::_rat_to_int;
		using base::_to_single_int;
		using base::_to_single_rat;
		using base::base;
		using base::err_tolerance;
		using base::log;
		using base::rgen;
		using base::threads;
		using base::tolerance;

		using Strategies = CutsGenStrategies<int_type, rational_type>;

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

	private:
		Hilbert_base_alg<intT, rationalT, lintT, lrationalT> _Hilbert_base_alg;

		void init_Hilbert_alg_params()
		{
			_Hilbert_base_alg.tolerance = tolerance;
			_Hilbert_base_alg.err_tolerance = err_tolerance;
			_Hilbert_base_alg.rgen = rgen;
			_Hilbert_base_alg.threads = threads;
		}

		is_mat _AT_copy;
		rs_vec _optLP_copy;

		std::vector<cut_info> _resulting_cuts;

		void copy_original_problem_info(const is_mat &A, const rs_vec &optLP)
		{
			_AT_copy = A;
			_optLP_copy = optLP;
		}

	protected:
		using transform_base = ILPCuts::transform_base<int_type, rational_type, lint_type, lrational_type>;
		using sparse_transform = ILPCuts::sparse_transform<int_type, rational_type, lint_type, lrational_type>;
		using lsparse_transform = ILPCuts::lsparse_transform<int_type, rational_type, lint_type, lrational_type>;
		struct dir_coeffs
		{
			rs_vec coeffs;
			size_t dir_ID;
		};

		struct cone_traits
		{
			is_mat A;
			lint_type det;
			std::vector<dir_coeffs> dir_coeffs_pull;

			void move_from_base_ptr(std::unique_ptr<cone_traits> base)
			{
				//*this = std::move(*base);
				A = std::move(base->A);
				det = std::move(base->det);
				dir_coeffs_pull = std::move(base->dir_coeffs_pull);
			}
		};

		virtual std::unique_ptr<cone_traits> compute_cone_traits(const cone_traits *C,
																 const transform_base &trans,
																 std::vector<dir_coeffs> dir_coeffs_pull)
		{

			size_t n = C->A.cols();
			size_t i = trans.get_i();

			auto new_cone = std::make_unique<cone_traits>();

			new_cone->dir_coeffs_pull = std::move(dir_coeffs_pull);

			try
			{
				out() << "Updating " << i << "-th column of A ...\n";
				new_cone->A = C->A;

				new_cone->A.col(i) = trans(C->A);

				Utility::prune_zeroes_int(new_cone->A);
			}
			catch (ToleranceErrorException exc)
			{
				out() << "A-column construction tolerance error: " << exc.what() << endl;
			}

			try
			{
				out() << "Updating a determinant ...\n";
				new_cone->det = trans[C->det];
			}
			catch (ToleranceErrorException exc)
			{
				out() << "Determinant recompute tolerance error: " << exc.what() << endl;
			}

			return std::move(new_cone);
		}

		std::vector<is_mat> partial_uni_decomp(std::unique_ptr<const cone_traits> C);

		std::vector<std::unique_ptr<cone_traits>> compute_layer(const cone_traits *C,
																std::vector<Atomic<size_t>> &directions_evaluation_counter);

		std::vector<dir_coeffs> compute_dirr_coeffs_pull(const cone_traits &C,
														 const transform_base &trans,
														 const std::vector<Atomic<size_t>> &directions_counter)
		{

			size_t n = C.A.cols();

			std::vector<dir_coeffs> res;

			for (size_t i = 0; i < C.dir_coeffs_pull.size(); ++i)
			{
				size_t dir_ID = C.dir_coeffs_pull[i].dir_ID;
				if (directions_counter[dir_ID] < cone_limit_per_direction)
				{

					const rs_vec &old_coeffs = C.dir_coeffs_pull[i].coeffs;
					rs_vec new_coeffs = trans.compute_coeffs(old_coeffs);

					if (new_coeffs.size() == n)
					{
						res.emplace_back(dir_coeffs{
							std::move(new_coeffs), dir_ID});
					}
				}
			}

			return std::move(res);
		}

		using process_func = std::function<void(const transform_base &)>;

		virtual void process_node(const cone_traits &C, process_func pf) = 0;

		void generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP);

		std::vector<rs_vec> generate_directions_pool(const is_mat &A)
		{
			std::vector<rs_vec> directions_pool;

			for (const auto &strategy : direction_strategies)
			{
				std::vector<rs_vec> strategy_directions = std::move(strategy->generate_directons(A));
				directions_pool.insert(directions_pool.end(),
									   std::make_move_iterator(strategy_directions.begin()),
									   std::make_move_iterator(strategy_directions.end()));
			}

			return std::move(directions_pool);
		}

		vector<cut_info> generate_cuts_from_row_matrix(const is_mat &cuts_mat, const rs_vec &optLP)
		{
			vector<cut_info> res;

			rs_vec rat_b = cuts_mat.template cast<rational_type>() * optLP;
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

	public:
		size_t cone_limit_per_direction = 1;
		std::vector<std::unique_ptr<typename Strategies::DirectionsStrategyBase>> direction_strategies;

		using size_func = std::function<size_t(size_t)>;
		size_func determinant_cut_limit = [](size_t n) -> size_t
		{ return n; };

		explicit ConicCutsGenerator_Base(logger &log)
			: base(log), _Hilbert_base_alg(log)
		{
		}

		logger::stream &out()
		{
			log.out() << "ILPCuts-Generator: ";
			return log.out();
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

		virtual void generate_cuts(const is_mat &A, const rs_vec &optLP) = 0;

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
			_AT_copy.resize(0, 0);
			_optLP_copy.resize(0);
		}

		void check_cuts_validity()
		{
			size_t n = _AT_copy.cols();
			assert(n > 0);

			for (size_t i = 0; i < _resulting_cuts.size(); ++i)
			{
				const auto &cut = _resulting_cuts[i];

				rs_vec t = Utility::system_solve(_AT_copy.template cast<rational_type>().eval(), cut.ray.template cast<rational_type>().eval());

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

				if (abs(cut.ray.template cast<rational_type>().dot(_optLP_copy) - cut.frac_rhs) > err_tolerance)
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

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::is_mat>
	ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::partial_uni_decomp(std::unique_ptr<const cone_traits> C)
	{
		size_t n = C->A.cols();
		omp_set_dynamic(std::max(threads, size_t(1)));
		std::vector<std::unique_ptr<const cone_traits>> current_layer;
		std::vector<std::vector<is_mat>> result_by_treads(threads);

		size_t directions_pool_size = C->dir_coeffs_pull.size();
		current_layer.push_back(std::move(C));

		while (!current_layer.empty())
		{
			std::vector<Atomic<size_t>> directions_evaluation_counter(directions_pool_size);

			size_t dynamic_threads = std::min(threads, current_layer.size());

			std::vector<std::vector<std::unique_ptr<const cone_traits>>> next_layers_by_treads(dynamic_threads);
			omp_set_num_threads(dynamic_threads);
#ifndef ILPCUTS_PARALLEL_OFF
#pragma omp parallel for
#endif
			for (size_t i = 0; i < current_layer.size(); ++i)
			{
				auto node = current_layer[i].get();
				out() << "Computing next layer ..." << endl;
				auto next_node_layer = compute_layer(node, directions_evaluation_counter);
				int tid = omp_get_thread_num();
				for (auto &next_node : next_node_layer)
				{
					if (next_node->det <= determinant_cut_limit(n))
					{
						out() << "Determinant is cutted on the value of " << determinant_cut_limit(n) << ",\n"
							  << "Total Gomory's cuts procedure is started..." << endl;
						result_by_treads[tid].push_back(_Hilbert_base_alg.compute(next_node->A, next_node->det));
					}
					else
						next_layers_by_treads[tid].push_back(std::move(next_node));
				}
			}

			std::swap(current_layer, next_layers_by_treads[0]);
			for (int i = 1; i < next_layers_by_treads.size(); ++i)
			{
				current_layer.insert(
					current_layer.end(),
					std::make_move_iterator(next_layers_by_treads[i].begin()),
					std::make_move_iterator(next_layers_by_treads[i].end()));
			}
		}

		std::vector<is_mat> &concatenated_result = result_by_treads[0];
		for (int i = 1; i < result_by_treads.size(); ++i)
		{
			concatenated_result.insert(
				concatenated_result.end(),
				std::make_move_iterator(result_by_treads[i].begin()),
				std::make_move_iterator(result_by_treads[i].end()));
		}

		return std::move(concatenated_result);
	}

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<std::unique_ptr<typename ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::cone_traits>>
	ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::compute_layer(const cone_traits *C,
																		  std::vector<Atomic<size_t>> &directions_evaluation_counter)
	{

		std::vector<std::unique_ptr<cone_traits>> result;

		const auto &A = C->A;
		const auto &det = C->det;
		size_t n = A.cols();

		out() << "Current determinant: " << det << "\n"
			  << "Nonzeros in A: " << A.nonZeros() << "\n"
			  << "A's inifinity norm: " << Utility::max_norm(A) << endl;

		auto process_new_candidate = [&](const transform_base &trans)
		{
			// out() << "Processing new cone candidate ...\n"
			// 	  << "Computing new directions coeff. pull ... \n";
			auto dirr_coeffs_pull = std::move(compute_dirr_coeffs_pull(
				*C, trans, directions_evaluation_counter));

			if (!dirr_coeffs_pull.empty())
			{
				// out() << "New candidate cone is approved ...\n";
				++directions_evaluation_counter[dirr_coeffs_pull[0].dir_ID];

				std::unique_ptr<cone_traits> new_cone = std::move(compute_cone_traits(
					C, trans, std::move(dirr_coeffs_pull)));

				result.push_back(std::move(new_cone));
			}
		};

		process_node(*C, process_new_candidate);

		return std::move(result);
	}

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	void ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP)
	{

		init_Hilbert_alg_params();

		clear_cuts();

		copy_original_problem_info(C->A, optLP);

		Timer timer;
		timer.start();

		out() << "Cone decomposition procedure is started ..." << endl;

		std::vector<is_mat> uni_mat_pool = std::move(partial_uni_decomp(
			std::move(C)));

		out() << "Cone decomposition procedure is finished, time = " << timer.finish_str();

		for (is_mat &uni_mat : uni_mat_pool)
		{
			uni_mat = uni_mat.transpose().eval();
		}

		out() << "Transforming cones to cuts..." << endl;
		timer.start();

		for (const is_mat &uni_mat : uni_mat_pool)
		{
			auto new_cuts = generate_cuts_from_row_matrix(uni_mat, optLP);

			_resulting_cuts.insert(_resulting_cuts.end(), make_move_iterator(new_cuts.begin()),
								   make_move_iterator(new_cuts.end()));
		}

		out() << "Transforming is finished, time = " << timer.finish_str()
			  << "Filtering duplicate cuts..." << endl;

		std::sort(_resulting_cuts.begin(), _resulting_cuts.end());
		auto last = std::unique(_resulting_cuts.begin(), _resulting_cuts.end());
		_resulting_cuts.erase(last, _resulting_cuts.end());
	}
}