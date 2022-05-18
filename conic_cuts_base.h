#pragma once

#include "base_module.h"

namespace ILPCuts
{

	template <typename intT, typename ratT>
	struct CutsGenStrategies
	{
		using int_type = intT;
		using rational_type = ratT;
		using is_mat = SparseMatrix<int_type>;
		using is_vec = SparseVector<int_type>;
		using rs_vec = SparseVector<rational_type>;

		class DirectionsStrategyBase
		{
		public:
			virtual std::vector<rs_vec> generate_directons(const is_mat &A) = 0;
		};

		class GreadyDirection : public DirectionsStrategyBase
		{
		private:
			rs_vec _c;

		public:
			GreadyDirection(const rs_vec &c)
				: _c(c) {}

			std::vector<rs_vec> generate_directons(const is_mat &A) override
			{
				std::vector<rs_vec> res(1);
				res[0] = _c;
				return std::move(res);
			}
		};

		class SomeDirection : public DirectionsStrategyBase
		{
		public:
			std::vector<rs_vec> generate_directons(const is_mat &A) override
			{
				size_t n = A.cols();
				rs_vec direction(n);
				direction.setZero();

				for (size_t i = 0; i < n; ++i)
				{
					direction += A.template cast<rational_type>().col(i);
				}
				direction /= rational_type(n);

				std::vector<rs_vec> res(1);
				res[0] = std::move(direction);
				return std::move(res);
			}
		};

		template <typename randGenT>
		class ConvUniform : public DirectionsStrategyBase
		{
		public:
			using gen_type = randGenT;

		private:
			const size_t _rand_dir_num;
			gen_type &_gen;

		public:
			ConvUniform(gen_type &gen, size_t rand_dir_num = 3)
				: DirectionsStrategyBase(), _gen(gen), _rand_dir_num(rand_dir_num) {}

			std::vector<rs_vec> generate_directons(const is_mat &A) override
			{
				std::vector<rs_vec> result_directions;
				size_t n = A.cols();

				result_directions.reserve(_rand_dir_num);

				std::uniform_real_distribution<double> udistr(0.0, 1.0);

				for (size_t i = 0; i < _rand_dir_num; ++i)
				{
					rs_vec new_dir = A.col(0).template cast<rational_type>();
					for (size_t j = 1; j < n; ++j)
					{
						rational_type alpha = udistr(_gen);
						new_dir = (new_dir * alpha + A.col(j).template cast<rational_type>() * (1.0 - alpha)).eval();
					}
					result_directions.push_back(std::move(new_dir));
				}

				return std::move(result_directions);
			}
		};
	};

	template <typename intT, typename rationalT, typename lintT, typename lrationalT>
	class ConicCutsGenerator_Base
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

		using Strategies = CutsGenStrategies<int_type, rational_type>;

	private:
		using me = ConicCutsGenerator_Base<intT, rationalT, lintT, lrationalT>;

		std::ostream *_out;
		std::ostream *_err;

	protected:
		using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;
		using LRatToLInt_Unary = Utility::RatToInt_Unary<lint_type, lrational_type>;
		std::unique_ptr<RatToInt_Unary> _rat_to_int;
		std::unique_ptr<LRatToLInt_Unary> _lrat_to_lint;

		Utility::ToSingle<lint_type, int_type> _to_single_int;
		Utility::ToSingle<lrational_type, rational_type> _to_single_rat;

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
				*this = std::move(*base);
				/*A = std::move(base->A);
				det = std::move(base->det);
				dir_coeffs_pull = std::move(base->dir_coeffs_pull);*/
			}
		};

		struct transform_base
		{
			me *owner;

			transform_base(me *aowner)
				: owner(aowner) {}

			virtual is_vec operator()(const is_mat &A) const = 0;
			virtual lint_type operator[](const lint_type &det) const = 0;
			virtual rs_vec compute_coeffs(const rs_vec &old_coefs) const = 0;
			virtual size_t get_i() const = 0;
			virtual lis_mat compute_adj(lis_mat old_adjA) const = 0;
		};

		struct sparce_transform : public transform_base
		{
			Utility::sparce_transform<int_type> data;

			sparce_transform(me *aowner, size_t ai,
							 const is_vec &at, const int_type &adenum)
				: transform_base(aowner), data{ai, at, adenum}
			{
			}

			is_vec operator()(const is_mat &A) const override
			{
				rs_vec res = (A * data.t).template cast<rational_type>() / data.denum;

				return res.unaryExpr(*(transform_base::owner->_rat_to_int));
			}

			lint_type operator[](const lint_type &det) const override
			{
				const int_type &ti = data.t.coeff(data.i);
				lrational_type res(static_cast<lrational_type>(det * ti) / data.denum);
				// cout << det << endl;
				// cout << (lint_type)(res * 2) << endl;
				return (*transform_base::owner->_lrat_to_lint)(res);
			}

			rs_vec compute_coeffs(const rs_vec &old_coefs) const override
			{
				return Utility::compute_transformed_coeffs(old_coefs,
														   data,
														   transform_base::owner->tolerance);
			}

			size_t get_i() const override { return data.i; }

			lis_mat compute_adj(lis_mat old_adjA) const override
			{
				size_t n = old_adjA.cols();
				is_mat Q = compute_adj_Q(n, data);
				old_adjA = ((Q.template cast<lint_type>() * old_adjA).template cast<lrational_type>() / data.denum)
							   .unaryExpr(*transform_base::owner->_lrat_to_lint);

				return std::move(old_adjA);
			}
		};

		struct lsparce_transform : public transform_base
		{
			Utility::sparce_transform<lint_type> data;

			lsparce_transform(me *aowner, size_t ai,
							  const lis_vec &at, const lint_type &adenum)
				: transform_base(aowner), data{ai, at, adenum}
			{
			}

			is_vec operator()(const is_mat &A) const override
			{
				lrs_vec res = (A.template cast<lint_type>() * data.t).template cast<lrational_type>() / data.denum;

				return res.unaryExpr(*(transform_base::owner->_lrat_to_lint)).unaryExpr(transform_base::owner->_to_single_int);
			}

			lint_type operator[](const lint_type &det) const override
			{
				const lint_type &ti = data.t.coeff(data.i);
				lrational_type res(static_cast<lrational_type>(det * ti) / data.denum);
				return (*transform_base::owner->_lrat_to_lint)(res);
			}

			rs_vec compute_coeffs(const rs_vec &old_coefs) const override
			{
				return Utility::compute_transformed_coeffs(old_coefs.template cast<lrational_type>().eval(),
														   data,
														   transform_base::owner->tolerance)
					.unaryExpr(transform_base::owner->_to_single_rat);
			}

			size_t get_i() const override { return data.i; }

			lis_mat compute_adj(lis_mat old_adjA) const override
			{
				size_t n = old_adjA.cols();
				lis_mat Q = compute_adj_Q(n, data);
				old_adjA = ((Q * old_adjA).template cast<lrational_type>() / data.denum)
							   .unaryExpr(*transform_base::owner->_lrat_to_lint);

				return std::move(old_adjA);
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
				new_cone->A = C->A;

				new_cone->A.col(i) = trans(C->A);

				Utility::prune_zeroes_int(new_cone->A);
			}
			catch (ToleranceErrorException exc)
			{
				out() << "A-column construction tolerance error: " << exc.what() << '\n';
				err() << "A-column construction tolerance error: " << exc.what() << '\n';
			}

			try
			{
				new_cone->det = trans[C->det];
			}
			catch (ToleranceErrorException exc)
			{
				out() << "Determinant recompute tolerance error: " << exc.what() << '\n';
				err() << "Determinant recompute tolerance error: " << exc.what() << '\n';
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

					const rs_vec &coeffs = C.dir_coeffs_pull[i].coeffs;

					rs_vec new_coeffs = trans.compute_coeffs(coeffs);

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

		std::vector<is_vec> generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP);

		std::vector<rs_vec> generate_directions_pool(const is_mat &A)
		{
			std::vector<rs_vec> directions_pool;

			for (auto &strategy : direction_strategies)
			{
				std::vector<rs_vec> strategy_directions = std::move(strategy->generate_directons(A));
				directions_pool.insert(directions_pool.end(), strategy_directions.begin(), strategy_directions.end());
			}

			return std::move(directions_pool);
		}

	public:
		double tolerance;
		double err_tolerance;
		size_t cone_limit_per_direction;
		size_t threads;
		std::vector<std::unique_ptr<typename Strategies::DirectionsStrategyBase>> direction_strategies;
		std::mt19937 *rgen;

		ConicCutsGenerator_Base(std::ostream *outs = &std::cout, std::ostream *errs = &std::cout)
			: _out(outs), _err(errs), cone_limit_per_direction(0),
			  tolerance(0), err_tolerance(0), threads(0), _rat_to_int(), rgen(nullptr)
		{
		}

		virtual void Init(double a_tolerance = 0.000001,
						  double a_err_tolerance = 0.001,
						  size_t a_cone_limit_per_direction = 1,
						  size_t a_threads = 8, std::mt19937 *a_rgen = nullptr)
		{

			tolerance = a_tolerance;
			err_tolerance = a_err_tolerance;
			cone_limit_per_direction = a_cone_limit_per_direction;
			// direction_strategies.push_back(std::make_unique<typename Strategies::SomeDirection>());
			threads = a_threads;
			rgen = a_rgen;

			_rat_to_int = std::make_unique<RatToInt_Unary>(err_tolerance);
			_lrat_to_lint = std::make_unique<LRatToLInt_Unary>(err_tolerance);
		}

		std::ostream &out() { return *_out; }

		std::ostream &err() { return *_err; }

		is_mat adjust_normal_cone(const is_mat &A, const is_vec &c)
		{
			is_mat resA = A;

			rs_vec t = Utility::system_solve(A.template cast<rational_type>().eval(), c.template cast<rational_type>().eval());
			/*cout << c << endl << endl;
			cout << t << endl;
			cin.ignore();*/

			for (typename rs_vec::InnerIterator it(t); it; ++it)
			{
				// it.value() < tolerance previouse code. working fast but its must be not correct
				if (it.value() < 0 && abs(it.value()) > tolerance)
				{
					/*cout << it.value() << endl;
					cin.ignore();*/
					resA.col(it.index()) = -resA.col(it.index());
				}
			}

			return std::move(resA);
		}

		virtual std::vector<is_vec> generate_cuts(const is_mat &A, const rs_vec &optLP) = 0;
	};

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::is_mat>
	ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::partial_uni_decomp(std::unique_ptr<const cone_traits> C)
	{
		omp_set_dynamic(std::max(threads, size_t(1)));
		std::vector<std::unique_ptr<const cone_traits>> current_layer;
		std::vector<std::vector<is_mat>> result(threads, std::vector<is_mat>());

		size_t directions_pool_size = C->dir_coeffs_pull.size();
		current_layer.push_back(std::move(C));

		while (!current_layer.empty())
		{

			/*out() << "Curr level size: " << current_layer.size() << std::endl;
			std::cin.ignore();*/

			std::vector<Atomic<size_t>> directions_evaluation_counter(directions_pool_size /*, Atomic<size_t>(0)*/);

			int dynamic_threads = std::min(threads, current_layer.size());

			std::vector<std::vector<std::unique_ptr<const cone_traits>>> next_layers(dynamic_threads /*, std::vector<cone_traits>()*/);
			omp_set_num_threads(dynamic_threads);
#ifndef ILPCUTS_PARALLEL_OFF
#pragma omp parallel for
#endif
			for (int i = 0; i < current_layer.size(); ++i)
			{
				auto node = current_layer[i].get();
				auto next_node_layer = std::move(compute_layer(node, directions_evaluation_counter));
				int tid = omp_get_thread_num();
				for (auto &next_node : next_node_layer)
				{
					if (next_node->det == 1)
						result[tid].push_back(std::move(next_node->A));
					else
						next_layers[tid].push_back(std::move(next_node));
				}
			}

			std::swap(current_layer, next_layers[0]);
			for (int i = 1; i < next_layers.size(); ++i)
			{
				current_layer.insert(
					current_layer.end(),
					std::make_move_iterator(next_layers[i].begin()),
					std::make_move_iterator(next_layers[i].end()));
			}
		}

		std::vector<is_mat> &concatenated_result = result[0];
		for (int i = 1; i < result.size(); ++i)
		{
			concatenated_result.insert(
				concatenated_result.end(),
				std::make_move_iterator(result[i].begin()),
				std::make_move_iterator(result[i].end()));
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

		out() << "Current determinant: " << det << std::endl;
		out() << "Nonzeros: " << A.nonZeros() << std::endl;
		out() << "A's inifinity norm: " << Utility::max_norm(A) << std::endl;

		auto process_new_candidate = [&](const transform_base &trans)
		{
			auto dirr_coeffs_pull = std::move(compute_dirr_coeffs_pull(
				*C, trans, directions_evaluation_counter));

			if (!dirr_coeffs_pull.empty())
			{
				++directions_evaluation_counter[dirr_coeffs_pull[0].dir_ID];

				std::unique_ptr<cone_traits> new_cone = std::move(compute_cone_traits(
					C, trans, dirr_coeffs_pull));

				result.push_back(std::move(new_cone));
				out() << "Adding step done..." << std::endl;
			}
		};

		process_node(*C, process_new_candidate);

		return std::move(result);
	}

	template <typename intT, typename rationalT, typename lintT, typename lratT>
	std::vector<typename ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::is_vec>
	ConicCutsGenerator_Base<intT, rationalT, lintT, lratT>::generate_cuts(std::unique_ptr<const cone_traits> C, const rs_vec &optLP)
	{

		auto start = std::chrono::steady_clock::now();

		std::vector<is_vec> res_cuts;

		std::vector<is_mat> uni_mat_pool = std::move(partial_uni_decomp(
			std::move(C)));

		for (is_mat &uni_mat : uni_mat_pool)
		{
			uni_mat = uni_mat.transpose();
		}

		for (const is_mat &uni_mat : uni_mat_pool)
		{
			rs_vec rat_b = uni_mat.template cast<rational_type>() * optLP;
			is_vec rounded_b = rat_b.unaryExpr(Utility::Round<rational_type>()).template cast<int_type>();
			is_vec floor_b = rat_b.unaryExpr(Utility::Floor<rational_type>()).template cast<int_type>();

			rs_vec diff_b = (rat_b - rounded_b.template cast<rational_type>()).cwiseAbs();

			/*std::cout << uni_mat << std::endl;
			std::cout << optLP.transpose() << std::endl;
			std::cout << rat_b.transpose() << std::endl;
			std::cout << rounded_b.transpose() << std::endl;
			std::cout << floor_b.transpose() << std::endl;
			std::cout << diff_b.transpose() << std::endl;*/

			typename rs_vec::InnerIterator diff_it(diff_b);
			// rs_vec::InnerIterator floor_it(floor_b);
			// size_t jjjj = 0;
			while (diff_it)
			{
				/*if (diff_it.value() != 0)
				{
					cout << ++jjjj << ", " << rat_b.coeff(diff_it.index()) - floor_b.coeff(diff_it.index()) << endl;
					cin.ignore();
				}*/
				if (diff_it.value() > rational_type(tolerance))
				{
					is_vec cut = uni_mat.row(diff_it.index());
					// cut.reserve(1);
					size_t n = cut.size();
					cut.conservativeResize(n + 1);
					cut.insert(n) = floor_b.coeff(diff_it.index());
					// cout << cut << endl;
					// out() << "cut has been added: (" << cut.transpose() << ")\n";
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
		res_cuts.erase(last, res_cuts.end());

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		out() << "ConicCuts: cuts generation time: " << elapsed_seconds.count() << "s\n";
		out() << "ConicCuts: total cuts: " << res_cuts.size() << "\n";

		return std::move(res_cuts);
	}
}