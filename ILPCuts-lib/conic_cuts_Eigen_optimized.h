#pragma once

#include"base_module.h"

namespace ILPCuts {

	using namespace Eigen;

	template<typename intT, typename ratT, typename lintT, typename lratT>
	class ConicCutsGeneratorOpt :public ConicCutsGenerator_Base<intT, ratT> {
	private:
		using base = typename ConicCutsGenerator_Base<intT, ratT>;
		using my_type = typename ConicCutsGeneratorOpt<intT, ratT, lintT, lratT>;

	public:
		using typename base::rational_type;
		using typename base::int_type;

		using typename base::is_mat;
		using typename base::rs_mat;

		using typename base::is_vec;
		using typename base::rs_vec;

		using lint_type = lintT;
		using lrational_type = lratT;

		using lis_mat = SparseMatrix<lint_type>;
		using lrs_mat = SparseMatrix<lrational_type>;

		using lis_vec = SparseVector<lint_type>;
		using lrs_vec = SparseVector<lrational_type>;

		using base::tolerance;
		using base::threads;
		using base::out;
		using base::err;
		using base::err_tolerance;
		using base::cone_limit_per_direction;
		using base::direction_strategies;
		using base::_rat_to_int;
		using base::rgen;

		using base::base;

	private:

		struct cone_traits {
			using decomp_type = SparseLU<rs_mat, COLAMDOrdering<int>>;
			is_mat A;
			lint_type det;
			std::unique_ptr<decomp_type> decomp;
			/*cone_traits(cone_traits&& a_cone) = default;
			cone_traits() = default;
			cone_traits(const cone_traits& a_cone) = default;*/
			
			//cone_traits& operator=(cone_traits&&) = default;
			//cone_traits(const cone_traits&) = delete;
			//cone_traits& operator=(const cone_traits&) = delete;
		};

		is_vec even_system_solve(const cone_traits& C) {
			size_t n = C.A.cols();

			is_vec res(n);

			bool odd_component_generated = false;
			while (!odd_component_generated) {

				//generation of a random rhs
				rs_vec rhs = Utility::bernouly_vec_generate<int_type>(n, *rgen);

				/*cout << rhs << endl;
				cin.ignore();*/

				//solve the system
				rs_vec x = std::move(Utility::system_solve(rhs, *C.decomp));

				/*cout << x << endl;
				cin.ignore();*/

				//create long integer solution
				lis_vec int_x = (x.cast<lrational_type>() * C.det).unaryExpr(Utility::RatToInt_Unary<lint_type, lrational_type>(err_tolerance));

				lint_type max_2_pow = Utility::vec_max2pow(int_x);
				if ((static_cast<lint_type>(C.det) / max_2_pow) % 2 == 0) {
					int_x /= max_2_pow;
					//out() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					odd_component_generated = true;

					/*cout << int_x << endl;
					cin.ignore();*/

					//find odd solutions and convert them to 1 in result vector
					res = int_x.unaryExpr(std::ref(Utility::mod2_unaryop<lint_type>)).cast<int_type>();

					/*cout << res << endl;
					cin.ignore();*/

					Utility::prune_zeroes_int(res);
					/*rs_vec Z = (A * res) / 2;
					cout << Z << endl;
					cin.ignore();*/

					/*cout << res << endl;
					cin.ignore();*/
				}
			}

			return std::move(res);
		}

		lis_vec odd_system_solve(const cone_traits& C) {
			size_t n = C.A.cols();

			lis_vec res(n);

			bool nonzero_component_generated = false;
			while (!nonzero_component_generated) {

				//generation of a random rhs
				rs_vec rhs = Utility::bernouly_vec_generate<int_type>(n, *rgen);

				//solve the system
				rs_vec x = std::move(Utility::system_solve(rhs, *C.decomp));

				//create long integer solution
				res = (x.cast<lrational_type>() * C.det).unaryExpr(Utility::RatToInt_Unary<lint_type, lrational_type>(err_tolerance)).cast<lint_type>();

				//transform modulo det
				res = res.unaryExpr(Utility::mod_unaryop<lint_type>(C.det));
				Utility::prune_zeroes_int(res);

				if (res.nonZeros() > 0) nonzero_component_generated = true;
			}

			return std::move(res);
		}

		template<typename funT>
		void process_even_case(const cone_traits& C, funT process_fun) {
			lis_vec x = even_system_solve(C).cast<lint_type>();
			for (typename lis_vec::InnerIterator it(x); it; ++it)
				if (it.value() > 0)
					process_fun(it.index(), x, 2);
		}

		template<typename funT>
		void process_odd_case(const cone_traits& C, funT process_fun) {
			lis_vec x = odd_system_solve(C);
			for (typename lis_vec::InnerIterator it(x); it; ++it)
				if (it.value() > 0) {
					lint_type lambda = Utility::mod_equation_solve<lint_type>(it.value(), C.det - gcd(C.det, it.value()), C.det);
					lis_vec t = (x * lambda).unaryExpr(Utility::mod_unaryop(C.det));

					process_fun(it.index(), t, C.det);
				}
		}

		class ProcessFunctor {
		private:
			my_type* _this;
			std::vector<cone_traits>& _res_cones;
			std::vector<Atomic<size_t>>& _directions_evaluation_counter;
			const cone_traits& _cone;
			const std::vector<rs_vec>& _directions_pool;

		public:
			ProcessFunctor(my_type* a_this, const cone_traits& a_cone, const std::vector<rs_vec>& a_directions_pool,
				std::vector<cone_traits>& a_res_cones, std::vector<Atomic<size_t>>& a_directions_evaluation_counter)
				:_this(a_this), _cone(a_cone), _directions_pool(a_directions_pool),
				_res_cones(a_res_cones), _directions_evaluation_counter(a_directions_evaluation_counter)
			{}

			void operator()(size_t idx, const lis_vec& t, const lint_type& denum) {

				const is_mat& A = _cone.A;
				const lint_type& det = _cone.det;
				const auto & decomp = _cone.decomp;
				size_t n = A.cols();

				cone_traits new_cone;

				try {
					/*if (t.coeff(i) == 0)
						std::cin.ignore();*/

					new_cone.A = A;
					new_cone.A.col(idx) = ((A.cast<lint_type>() * t).cast<lrational_type>() / denum)
						.unaryExpr(Utility::RatToInt_Unary<int_type, lrational_type>(_this->err_tolerance));

					Utility::prune_zeroes_int(new_cone.A);

					new_cone.decomp = std::make_unique<cone_traits::decomp_type>(new_cone.A.cast<rational_type>());

					//std::cout << A << std::endl;
					//std::cout << new_cone.A << std::endl;
					//.round()
					//.cast<int_type>();
				}
				catch (ToleranceErrorException exc) {
					_this->out() << "b-column construction tolerance error: " << exc.what() << '\n';
					_this->err() << "b-column construction tolerance error: " << exc.what() << '\n';
				}

				try {

					lint_type ti = t.coeff(idx);
					new_cone.det = Utility::rat_to_int<lint_type, lrational_type>
						(static_cast<lrational_type>(det * ti) / denum, _this->err_tolerance );

					/*if (new_cone.det == 0) {
						out() << det <<", "<< ti << std::endl;
						std::cin.ignore();
					}*/

				}
				catch (ToleranceErrorException exc) {
					_this->out() << "Determinant/inverse computing tolerance error: " << exc.what() << '\n';
					_this->err() << "Determinant/inverse computing tolerance error: " << exc.what() << '\n';
				}

				for (size_t i = 0; i < _directions_pool.size(); ++i) {

					if (_directions_evaluation_counter[i] < _this->cone_limit_per_direction) {

						rs_vec duals = Utility::system_solve(_directions_pool[i], *new_cone.decomp);
						//rat_vec duals = adjA.cast<rational_type>() * directions_pool.col(i);
						/*cout << "duals: " << duals << endl;
						cin.ignore();*/

						/*int_mat I = A * adjA;
						cout << I << endl;
						cin.ignore();*/

						bool is_nontrivial_nonnegative = false;
						for (typename rs_vec::InnerIterator it(duals); it; ++it) {
							if (it.value() + _this->tolerance < 0) {
								is_nontrivial_nonnegative = false;
								break;
							}
							else if (it.value() + _this->tolerance > 0) {
								is_nontrivial_nonnegative = true;
							}
						}

						if (is_nontrivial_nonnegative) {
							if (_directions_evaluation_counter[i] >= _this->cone_limit_per_direction)
								break;
							++_directions_evaluation_counter[i];
							_res_cones.push_back(std::move(new_cone));
							_this->out() << "adding step done" << std::endl;
							break;
						}
					}

					//if (((duals + rational_type(tolerance)) > 0).all()) {
					//	result.push_back(std::move(new_cone));
					//	out() << "adding step done" << std::endl;
					//	//cin.ignore();
					//}
				}
			}
		};

		friend class ProcessFunctor;

		std::vector<cone_traits> process_cone(const cone_traits& cone, const std::vector<rs_vec>& directions_pool,
			std::vector<Atomic<size_t>>& directions_evaluation_counter) {
			std::vector<cone_traits> res_cones;

			const is_mat& A = cone.A;
			const lint_type& det = cone.det;
			size_t n = A.cols();

			out() << "Current determinant: " << det << std::endl;
			out() << "Nonzeros: " << A.nonZeros() << std::endl;

			ProcessFunctor process_functor(this, cone, directions_pool, res_cones, directions_evaluation_counter);
			if (det % 2 == 0) {
				out() << "even step..." << endl;
				process_even_case(cone, process_functor);
			}
			else {
				out() << "odd step..." << endl;
				process_odd_case(cone, process_functor);
			}

			/*cout << res_cones.size() << endl;
			cin.ignore();*/

			return std::move(res_cones);
		}

		std::vector<is_mat> partial_uni_decomp
		(const is_mat& A, const std::vector<rs_vec>& directions_pool,
			const lint_type& det) {

			omp_set_dynamic(max(threads, size_t(1)));
			std::vector<cone_traits> current_layer;
			std::vector<std::vector<is_mat>> result(threads, std::vector<is_mat>());

			current_layer.emplace_back(cone_traits{ A, det, 
				std::make_unique<cone_traits::decomp_type>(A.cast<rational_type>()) });

			while (!current_layer.empty()) {

				/*out() << "Curr level size: " << current_layer.size() << std::endl;
				std::cin.ignore();*/

				std::vector<Atomic<size_t>> directions_evaluation_counter(directions_pool.size(), Atomic<size_t>(0));

				int dynamic_threads = std::min(threads, current_layer.size());

				std::vector<std::vector<cone_traits>> next_layers(dynamic_threads);
				omp_set_num_threads(dynamic_threads);

#pragma omp parallel for
				for (int i = 0; i < current_layer.size(); ++i) {
					const cone_traits& node = current_layer[i];
					std::vector<cone_traits> next_node_layer = std::move(process_cone(node, directions_pool, directions_evaluation_counter));
					int tid = omp_get_thread_num();
					for (cone_traits& next_node : next_node_layer) {
						if (next_node.det == 1) {
							result[tid].push_back(std::move(next_node.A));
						} else {
							next_layers[tid].push_back(std::move(next_node));
						}
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

	public:

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP, const lint_type& det) {

			auto start = std::chrono::steady_clock::now();

			std::vector<is_vec> res_cuts;

			std::vector<rs_vec> directions_pool;
			for (auto& strategy : direction_strategies) {
				std::vector<rs_vec> strategy_directions = std::move(strategy->generate_directons(A));
				//out() << strategy_directions[0] << std::endl;
				directions_pool.insert(directions_pool.end(), strategy_directions.begin(), strategy_directions.end());
			}

			//auto directions_pool = std::move(directions_strategy.generate_directons(A));

			std::vector<is_mat> uni_mat_pool = std::move(partial_uni_decomp(A, directions_pool, det));

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
				while (diff_it) {
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

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP) {
			lint_type det = Utility::compute_abs_det(A);

			return std::move(generate_cuts(A, optLP, det));
		}
	};
}
