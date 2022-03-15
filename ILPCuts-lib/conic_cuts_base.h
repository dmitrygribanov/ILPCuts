#pragma once

#include "base_module.h"

namespace ILPCuts {

	template<typename intT, typename ratT>
	struct CutsGenStrategies {
		using int_type = intT;
		using rational_type = ratT;
		using is_mat = SparseMatrix<int_type>;
		using is_vec = SparseVector<int_type>;
		using rs_vec = SparseVector<rational_type>;

		class DirectionsStrategyBase {
		public:
			virtual std::vector<rs_vec> generate_directons(const is_mat& A) = 0;
		};

		class GreadyDirection :public DirectionsStrategyBase {
		private:
			rs_vec _c;

		public:
			GreadyDirection(const rs_vec& c)
				:_c(c) {}

			std::vector<rs_vec> generate_directons(const is_mat& A) override {
				std::vector<rs_vec> res(1);
				res[0] = _c;
				return std::move(res);
			}
		};

		class SomeDirection :public DirectionsStrategyBase {
		public:

			std::vector<rs_vec> generate_directons(const is_mat& A) override {
				size_t n = A.cols();
				rs_vec direction(n);
				direction.setZero();

				for (size_t i = 0; i < n; ++i) {
					direction += A.cast<rational_type>().col(i);
				}
				direction /= n;

				std::vector<rs_vec> res(1);
				res[0] = std::move(direction);
				return std::move(res);
			}
		};

		template<typename randGenT>
		class ConvUniform :public DirectionsStrategyBase {
		public:
			using gen_type = randGenT;

		private:
			const size_t _rand_dir_num;
			gen_type& _gen;

		public:

			ConvUniform(gen_type& gen, size_t rand_dir_num = 3)
				:DirectionsStrategyBase(), _gen(gen), _rand_dir_num(rand_dir_num) {}

			std::vector<rs_vec> generate_directons(const is_mat& A) override {
				std::vector<rs_vec> result_directions;
				size_t n = A.cols();

				result_directions.reserve(_rand_dir_num);

				std::uniform_real_distribution<double> udistr(0.0, 1.0);

				for (size_t i = 0; i < _rand_dir_num; ++i) {
					rs_vec new_dir = A.col(0).cast<rational_type>();
					for (size_t j = 1; j < n; ++j) {
						rational_type alpha = udistr(_gen);
						new_dir = (new_dir * alpha + A.col(j).cast<rational_type>() * (1.0 - alpha)).eval();
					}
					result_directions.push_back(std::move(new_dir));
				}

				return std::move(result_directions);
			}
		};
	};

	template<typename intT, typename rationalT, typename lintT, typename lrationalT>
	class ConicCutsGenerator_Base {
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


		using Strategies = typename CutsGenStrategies<int_type, rational_type>;

	private:

		std::ostream* _out;
		std::ostream* _err;

	protected:
		using RatToInt_Unary = typename Utility::RatToInt_Unary<int_type, rational_type>;
		using LRatToLInt_Unary = typename Utility::RatToInt_Unary<lint_type, lrational_type>;
		std::unique_ptr<RatToInt_Unary> _rat_to_int;
		std::unique_ptr<LRatToLInt_Unary> _lrat_to_lint;

		Utility::ToSingle<lint_type, int_type> _to_single_int;
		Utility::ToSingle<lrational_type, rational_type> _to_single_rat;

		struct dir_coeffs {
			rs_vec coeffs;
			size_t dir_ID;
		};

		struct cone_traits {
			is_mat A;
			lint_type det;
			std::vector<dir_coeffs> dir_coeffs_pull;
		};

		using sparce_transform = Utility::sparce_transform<lint_type>;

		std::vector<dir_coeffs> compute_dirr_coeffs_pull(const cone_traits& C,
			const sparce_transform& trans,
			const std::vector<Atomic<size_t>>& directions_counter) {

			size_t n = C.A.cols();

			std::vector<dir_coeffs> res;

			for (size_t i = 0; i < C.dir_coeffs_pull.size(); ++i) {
				size_t dir_ID = C.dir_coeffs_pull[i].dir_ID;
				if (directions_counter[dir_ID] < cone_limit_per_direction) {

					const rs_vec& coeffs = C.dir_coeffs_pull[i].coeffs;

					rs_vec new_coeffs = std::move(Utility::compute_transformed_coeffs(
						coeffs,
						trans,
						tolerance
					));

					if (new_coeffs.size() == n) {
						res.emplace_back(dir_coeffs{
							std::move(new_coeffs), dir_ID
							});
					}
				}
			}

			return std::move(res);
		}

	public:

		double tolerance;
		double err_tolerance;
		size_t cone_limit_per_direction;
		size_t threads;
		std::vector<std::unique_ptr<typename Strategies::DirectionsStrategyBase> > direction_strategies;
		std::mt19937* rgen;

		ConicCutsGenerator_Base(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:_out(outs), _err(errs), cone_limit_per_direction(0),
			tolerance(0), err_tolerance(0), threads(0), _rat_to_int(), rgen(nullptr)
		{}

		virtual void Init(double a_tolerance = 0.000001,
			double a_err_tolerance = 0.001,
			size_t a_cone_limit_per_direction = 1,
			size_t a_threads = 8, std::mt19937* a_rgen = nullptr) {

			tolerance = a_tolerance;
			err_tolerance = a_err_tolerance;
			cone_limit_per_direction = a_cone_limit_per_direction;
			direction_strategies.push_back(std::make_unique<Strategies::SomeDirection>());
			threads = a_threads;
			rgen = a_rgen;

			_rat_to_int = std::make_unique<RatToInt_Unary>(err_tolerance);
			_lrat_to_lint = std::make_unique<LRatToLInt_Unary>(err_tolerance);
		}

		std::ostream& out() { return *_out; }

		std::ostream& err() { return *_err; }

		is_mat adjust_normal_cone(const is_mat& A, const is_vec& c) {
			is_mat resA = A;

			rs_vec t = Utility::system_solve(A.cast<rational_type>().eval(), c.cast<rational_type>().eval());
			/*cout << c << endl << endl;
			cout << t << endl;
			cin.ignore();*/

			for (typename rs_vec::InnerIterator it(t); it; ++it) {
				// it.value() < tolerance previouse code. working fast but its must be not correct
				if (it.value() < 0 && abs(it.value()) > tolerance) {
					/*cout << it.value() << endl;
					cin.ignore();*/
					resA.col(it.index()) = -resA.col(it.index());
				}
			}

			return std::move(resA);
		}

		virtual std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP) = 0;

	};

}