#pragma once

#include<Eigen/SparseCore>
#include<Eigen/Dense>
#include<Eigen/SparseLU>

#include <mpirxx.h>

#include<vector>
#include<numeric>
#include<functional>
#include<sstream>
#include<random>
#include<chrono>

#include"base_module.h"

namespace ILPCuts {

	using namespace Eigen;

	template<typename intT, typename rationalT>
	class ConicCutsGenerator {
	public:
		using rational_type = rationalT;
		using int_type = intT;

		using is_mat = SparseMatrix<int_type>;
		using rs_mat = SparseMatrix<rational_type>;

		using is_vec = SparseVector<int_type>;
		using rs_vec = SparseVector<rational_type>;

		using id_mat = Matrix<int_type, Dynamic, Dynamic>;
		using rd_mat = Matrix<rational_type, Dynamic, Dynamic>;

		using id_vec = Matrix<int_type, Dynamic, 1>;
		using rd_vec = Matrix<rational_type, Dynamic, 1>;

	private:

		std::ostream* _out;
		std::ostream* _err;

		static int_type mod2_unaryop(const int_type& x) {
			return abs(x % 2);
		}

		class mod_unaryop {
		private:
			const int_type& _mod;

		public:
			mod_unaryop(const int_type& mod)
				:_mod(mod) {}

			int_type operator()(const int_type& x) const {
				int_type res = x % _mod;
				if (res < 0) res += _mod;
				return res;
			}
		};

		static int_type rat_to_int(const rational_type& rat_x, const rational_type& err_tolerance);

		class Round {
		public:
			rational_type operator()(const rational_type& rat_x) const {
				return std::move(round(rat_x));
			}
		};

		class Floor {
		public:
			rational_type operator()(const rational_type& rat_x) const {
				return std::move(floor(rat_x));
			}
		};

		class RatToInt_Unary {
		private:
			rational_type _err_tolerance;

		public:
			RatToInt_Unary(const rational_type& err_tolerance)
				:_err_tolerance(err_tolerance) {}

			int_type operator()(const rational_type& rat_x) const {
				return rat_to_int(rat_x, _err_tolerance);
			}
		};

		class ToleranceErrorException : public std::exception {
		public:
			ToleranceErrorException()
				:std::exception() {}

			ToleranceErrorException(std::string msg)
				:std::exception(msg.c_str()) {}
		};

		template<typename solverT>
		rs_mat compute_inverse(const is_mat & A, solverT & solver);

		std::tuple<is_mat, int_type> compute_adjugate_and_absdet(const is_mat& A);

		std::vector<is_mat> partial_uni_decomp
		(const is_mat& A, const std::vector<rs_vec>& directions_pool,
			const is_mat& adjA, const int_type& det);

		struct cone_traits {
			is_mat A;
			is_mat adjA;
			int_type det;
		};

		std::vector<cone_traits> compute_layer(const cone_traits& cone, const std::vector<rs_vec>& directions_pool, std::vector<size_t>& directions_evaluation_counter);

		int_type vec_max2pow(const is_vec& vec) const;

		int_type _mod_equation_solve(int_type a, int_type b, int_type mod) const;

		template<typename vecT>
		class SparseVecComp {
		public:
			bool operator()(const vecT& v1, const vecT& v2) {
				typename vecT::InnerIterator it1(v1);
				typename vecT::InnerIterator it2(v2);
				while (it1 && it2) {
					if (it1.index() > it2.index()) {
						if (it1.value() < 0)
							return true;
						else
							return false;
					}
					else if (it1.index() < it2.index()) {
						if (it2.value() > 0)
							return true;
						else
							return false;
					}
					else {
						if (it1.value() < it2.value())
							return true;
						else if (it1.value() > it2.value())
							return false;
					}

					++it1; ++it2;
				}

				if (it1) {
					if (it1.value() < 0)
						return true;
					else
						return false;
				}

				if (it2) {
					if (it2.value() > 0)
						return true;
					else
						return false;
				}

				return false;
			}
		};

		template<typename vecT>
		class SparseVecEquiv {
		public:
			bool operator()(const vecT& v1, const vecT& v2) {
				typename vecT::InnerIterator it1(v1);
				typename vecT::InnerIterator it2(v2);

				while (it1 && it2) {
					
					if (it1.index() != it2.index())
						return false;
					else {
						if (it1.value() != it2.value())
							return false;
					}

					++it1; ++it2;
				}

				return true;
			}
		};

	public:

		class DirectionsStrategyBase {
		public:
			virtual std::vector<rs_vec> generate_directons(const is_mat& A) = 0; // { return rat_mat();  };
		};

		rational_type tolerance;
		rational_type err_tolerance;
		size_t cone_limit_per_direction;
		size_t threads;
		std::vector<std::unique_ptr<DirectionsStrategyBase> > direction_strategies;

		ConicCutsGenerator(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:_out(outs), _err(errs), cone_limit_per_direction(0), 
			tolerance(0), err_tolerance(0), threads(0) {}

		void Init(const rational_type & a_tolerance = 0.00001,
			const rational_type& a_err_tolerance = 0.001, 
			size_t a_cone_limit_per_direction = 1,
			size_t a_threads = 8) {
			tolerance = a_tolerance;
			err_tolerance = a_err_tolerance;
			cone_limit_per_direction = a_cone_limit_per_direction;
			direction_strategies.push_back(std::make_unique<SomeDirection>() );
			threads = a_threads;
		}

		std::ostream& out() { return *_out; }

		std::ostream& err() { return *_err; }

		class GreadyDirection :public DirectionsStrategyBase {
		private:
			rs_vec _c;

		public:
			GreadyDirection(const rs_vec& c)
				:_c(c) {}

			std::vector<rs_vec> generate_directons(const is_mat & A) {
				std::vector<rs_vec> res(1);
				res[0] = _c;
				return std::move(res);
			}
		};

		class SomeDirection :public DirectionsStrategyBase {
		private:
			//rs_vec _direction;

		public:
			/*SomeDirection(const is_mat& A) {
				size_t n = A.cols();
				_direction.resize(n);
				_direction.setZero();

				for (size_t i = 0; i < n; ++i) {
					_direction += A.cast<rational_type>().col(i);
				}
				_direction /= n;
			}*/

			std::vector<rs_vec> generate_directons(const is_mat& A) {
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
			//const is_mat& _A;
			const size_t _rand_dir_num;
			gen_type& _gen;

		public:

			/*ConvUniform(gen_type& gen, const is_mat& A, size_t rand_dir_num = 3)
				:DirectionsStrategyBase(), _gen(gen), _A(A), _rand_dir_num(rand_dir_num) {}*/

			ConvUniform(gen_type& gen, size_t rand_dir_num = 3)
				:DirectionsStrategyBase(), _gen(gen), _rand_dir_num(rand_dir_num) {}

			std::vector<rs_vec> generate_directons(const is_mat& A) {
				std::vector<rs_vec> result_directions;
				size_t n = A.cols();

				result_directions.reserve(_rand_dir_num);

				//rat_mat NA = _A.cast<rational_type>().colwise().normalize();
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

		std::vector<is_vec> generate_cuts(const is_mat& A, const is_mat& adjA, const rs_vec& optLP, const int_type& det);

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP);

		is_mat adjust_normal_cone(const is_mat& A, const rs_vec& c);
	};

	template<typename intT, typename ratT>
	typename ConicCutsGenerator<intT, ratT>::is_mat 
		ConicCutsGenerator<intT, ratT>::adjust_normal_cone(const is_mat& A, const rs_vec& c) {

		is_mat resA = A;

		is_mat adjA;
		int_type det;

		std::tie(adjA, det) = compute_adjugate_and_absdet(A);

		rs_vec t = adjA.cast<rational_type>() * c;

		for (typename rs_vec::InnerIterator it(t); it; ++it) {
			if (it.value()  < tolerance) {
				resA.col(it.index()) = -resA.col(it.index());
			}
		}

		/*out() << t.transpose() << std::endl;
		out() << A << std::endl;
		out() << "--------------------------------------" << std::endl;
		out() << resA << std::endl;

		std::tie(adjA, det) = compute_adjugate_and_absdet(resA);
		t = adjA.cast<rational_type>() * c;

		out() << t.transpose() << std::endl;
		out() << c.transpose() << std::endl;*/

		return std::move(resA);
	}

	template<typename intT, typename ratT>
	std::vector<typename ConicCutsGenerator<intT, ratT>::is_mat> 
		ConicCutsGenerator<intT, ratT>::partial_uni_decomp
		(const is_mat& A, const std::vector<rs_vec>& directions_pool,
			const is_mat& adjA, const int_type& det) {
		std::vector<cone_traits> current_layer;
		std::vector<cone_traits> next_layer;
		std::vector<is_mat> result;

		current_layer.push_back(std::move(cone_traits{ A,adjA,det }));

		while (!current_layer.empty()) {

			/*out() << "Curr level size: " << current_layer.size() << std::endl;
			std::cin.ignore();*/

			std::vector<size_t> directions_evaluation_counter(directions_pool.size(), 0);

			for (const cone_traits& node : current_layer) {
				std::vector<cone_traits> next_node_layer = std::move(compute_layer(node, directions_pool, directions_evaluation_counter));
				for (cone_traits& next_node : next_node_layer) {
					if (next_node.det == 1)
						result.push_back(std::move(next_node.A));
					else
						next_layer.push_back(std::move(next_node));
				}
			}

			std::swap(current_layer, next_layer);
			next_layer.clear();
		}

		return std::move(result);
	}

	template<typename intT, typename ratT>
	std::vector<typename ConicCutsGenerator<intT, ratT>::is_vec>
		ConicCutsGenerator<intT, ratT>::generate_cuts(const is_mat& A, const is_mat& adjA, const rs_vec& optLP, const int_type& det) {

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
			is_vec rounded_b = rat_b.unaryExpr(Round()).cast<int_type>();
			is_vec floor_b = rat_b.unaryExpr(Floor()).cast<int_type>();

			rs_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).cwiseAbs();

			/*std::cout << uni_mat << std::endl;
			std::cout << optLP.transpose() << std::endl;
			std::cout << rat_b.transpose() << std::endl;
			std::cout << rounded_b.transpose() << std::endl;
			std::cout << floor_b.transpose() << std::endl;
			std::cout << diff_b.transpose() << std::endl;*/

			typename rs_vec::InnerIterator diff_it(diff_b);
			//rs_vec::InnerIterator floor_it(floor_b);
			while(diff_it) {
				if (diff_it.value() > tolerance) {
					is_vec cut = uni_mat.row(diff_it.index());
					cut.reserve(1);
					cut.insert(cut.size() - 1) = floor_b.coeff(diff_it.index());
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

		std::sort(res_cuts.begin(), res_cuts.end(), SparseVecComp<rs_vec>());
		auto last = std::unique(res_cuts.begin(), res_cuts.end(), SparseVecEquiv<rs_vec>());
		// v now holds {1 2 3 4 5 x x}, where 'x' is indeterminate
		res_cuts.erase(last, res_cuts.end());

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		out() << "ConicCuts: cuts generation time: " << elapsed_seconds.count() << "s\n";
		out() << "ConicCuts: total cuts: " << res_cuts.size() << "\n";

		return std::move(res_cuts);
	}

	template<typename intT, typename ratT>
	std::vector<typename ConicCutsGenerator<intT, ratT>::is_vec>
		ConicCutsGenerator<intT, ratT>::generate_cuts(const is_mat& A, const rs_vec& optLP) {

		//SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU(A.cast<rational_type>());

		int_type det = 0;
		is_mat adjA;

		//std::cout << A << std::endl;

		try {
			//det = rat_to_int(abs(decomp_LU.determinant()), err_tolerance);
			auto start = std::chrono::steady_clock::now();
			std::tie(adjA, det) = std::move(compute_adjugate_and_absdet(A));
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			out() << "ConicCuts: matrix inversion time: " << elapsed_seconds.count() << "s\n";

			return std::move(generate_cuts(A, adjA, optLP, det));
		}
		catch (ToleranceErrorException exc) {
			out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
		}

		//try {
		//	is_mat adjA = compute_inverse();
		//	// = (decomp_LU.inverse().array() * det).unaryExpr(RatToInt_Unary(err_tolerance));

		//	return std::move(generate_cuts(A, adjA, optLP, directions_strategy,
		//		det, tolerance, err_tolerance));
		//}
		//catch (ToleranceErrorException exc) {
		//	out() << "Root: inverse tolerance error: " << exc.what() << '\n';
		//	err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
		//}
	}

	template<typename intT, typename ratT>
	typename ConicCutsGenerator<intT, ratT>::int_type ConicCutsGenerator<intT, ratT>::rat_to_int(const rational_type& rat_x, const rational_type& err_tolerance) {
		//std::cout << rat_x << std::endl;
		int_type int_x = static_cast<int_type>(round(rat_x));
		/*std::cout << int_x << std::endl;
		std::cout << abs(rat_x - int_x) << std::endl;*/
		if (abs(rat_x - int_x) <= err_tolerance) {
			return int_x;
		}
		else {
			std::stringstream msg;
			msg << "rat_to_int tolerance error (" << abs(rat_x - int_x) << ")";
			//std::cout << abs(rat_x - int_x) << std::endl;
			throw ToleranceErrorException(msg.str());
		}
		return 0;
	}

	template<typename intT, typename ratT>
	std::tuple<typename ConicCutsGenerator<intT, ratT>::is_mat, typename ConicCutsGenerator<intT, ratT>::int_type> 
		ConicCutsGenerator<intT, ratT>::compute_adjugate_and_absdet(const is_mat& A) {
		//SparseLU<rs_mat> decomp_LU(A.cast<rational_type>());
		SparseLU<rs_mat, COLAMDOrdering<int>> decomp_LU(A.cast<rational_type>());

		int_type det = rat_to_int(abs(decomp_LU.determinant()), err_tolerance);
		is_mat adjA = (compute_inverse(A, decomp_LU) * det).unaryExpr(RatToInt_Unary(err_tolerance));

		size_t n = A.cols();
		rs_mat I(n, n);
		I.reserve(VectorXi::Constant(n, 2));
		for (size_t i = 0; i < n; ++i)
			I.insert(i, i) = det;

		/*if (rs_mat(A * adjA) == I)
			out() << "aaaaaaaaaaaaaaaaaaaa!!!!" << std::endl;*/

		return std::make_tuple(std::move(adjA), std::move(det));
	}


	template<typename intT, typename ratT>
	template<typename solverT>
	typename ConicCutsGenerator<intT, ratT>::rs_mat 
		ConicCutsGenerator<intT, ratT>::compute_inverse(const is_mat& A, solverT & solver) {

		size_t n = A.cols();

		rs_mat I(n, n);
		I.reserve(VectorXi::Constant(n, 2));
		for (size_t i = 0; i < n; ++i)
			I.insert(i, i) = 1;

		rs_mat res(n, n);
		for (size_t i = 0; i < n; ++i) {
			res.col(i) = solver.solve(I.col(i));
		}

		return std::move(res);
	}

	template<typename intT, typename ratT>
	std::vector<typename ConicCutsGenerator<intT, ratT>::cone_traits> 
		ConicCutsGenerator<intT, ratT>::compute_layer(const cone_traits& cone, const std::vector<rs_vec>& directions_pool, 
			std::vector<size_t> & directions_evaluation_counter) {

		std::vector<cone_traits> result;

		const is_mat& A = cone.A;
		const is_mat& adjA = cone.adjA;
		const int_type& det = cone.det;
		size_t n = A.cols();

		out() << "current determinant: " << det << std::endl;
		out() << "Nonzeros: " << A.nonZeros() << std::endl;

		auto process_new_candidate = [&](size_t i, const is_vec& t, const int_type& div)
		{
			cone_traits new_cone;

			try {
				/*if (t.coeff(i) == 0)
					std::cin.ignore();*/

				new_cone.A = A;
				new_cone.A.col(i) = ((A * t).cast<rational_type>() / div)
					.unaryExpr(RatToInt_Unary(err_tolerance));

				//std::cout << A << std::endl;
				//std::cout << new_cone.A << std::endl;
				//.round()
				//.cast<int_type>();
			}
			catch (ToleranceErrorException exc) {
				out() << "b-column construction tolerance error: " << exc.what() << '\n';
				err() << "b-column construction tolerance error: " << exc.what() << '\n';
			}

			//PartialPivLU<rat_mat> decomp_LU(new_cone.A.cast<rational_type>());

			try {

				int_type ti = t.coeff(i);
				new_cone.det = rat_to_int(static_cast<rational_type>(det * ti) / div, err_tolerance);

				/*if (new_cone.det == 0) {
					out() << det <<", "<< ti << std::endl;
					std::cin.ignore();
				}*/

				is_mat MM(n,n);
				using Tripl = Triplet<int_type>;
				std::vector<Tripl> triplets;
				for (size_t j = 0; j < n; ++j) {
					if( j != i )
						triplets.push_back(Tripl(j, j, ti));
					else
						triplets.push_back(Tripl(j, j, div));
				}

				for (typename is_vec::InnerIterator it(t); it; ++it) {
					if( it.value() > 0 && it.index() != i)
						triplets.push_back(Tripl(it.index(), i, -it.value()));
				}

				MM.setFromTriplets(triplets.begin(), triplets.end() );

				new_cone.adjA = ((MM * adjA).cast<rational_type>() / div).unaryExpr(RatToInt_Unary(err_tolerance));


				//std::tie(new_cone.adjA, new_cone.det) = compute_adjugate_and_absdet(new_cone.A);
				/*if (det % 2 == 0) {
					cout << "old det:" << det << endl;
					cout << "new det:" << new_cone.det << endl;
				}*/
			}
			catch (ToleranceErrorException exc) {
				out() << "Determinant/inverse computing tolerance error: " << exc.what() << '\n';
				err() << "Determinant/inverse computing tolerance error: " << exc.what() << '\n';
			}

			/*try {
				new_cone.adjA = (decomp_LU.inverse().array() * new_cone.det)
					.round()
					.cast<int_type>();
			}
			catch (ToleranceErrorException exc) {
				out() << "Inversion tolerance error: " << exc.what() << '\n';
				err() << "Inversion tolerance error: " << exc.what() << '\n';
			}*/

			for (size_t i = 0; i < directions_pool.size(); ++i) {
				//cout << "direction: " << directions_pool.col(i).transpose() << endl;
				/*if (directions_evaluation_counter[i] >= 3)
				{
					out() << directions_evaluation_counter[i] << std::endl;
					std::cin.ignore();
				}*/
				if (directions_evaluation_counter[i] < cone_limit_per_direction) {

					rs_vec duals = new_cone.adjA.cast<rational_type>() * directions_pool[i];
					//rat_vec duals = adjA.cast<rational_type>() * directions_pool.col(i);
					//cout << "duals: " << duals << endl;

					/*int_mat I = A * adjA;
					cout << I << endl;
					cin.ignore();*/

					bool is_nontrivial_nonnegative = false;
					for (typename rs_vec::InnerIterator it(duals); it; ++it) {
						if (it.value() + tolerance < 0) {
							is_nontrivial_nonnegative = false;
							break;
						}
						else if (it.value() + tolerance > 0 ) {
							is_nontrivial_nonnegative = true;
						}
					}

					if (is_nontrivial_nonnegative) {
						++directions_evaluation_counter[i];
						result.push_back(std::move(new_cone));
						out() << "adding step done" << std::endl;
						break;
					}
				}

				//if (((duals + rational_type(tolerance)) > 0).all()) {
				//	result.push_back(std::move(new_cone));
				//	out() << "adding step done" << std::endl;
				//	//cin.ignore();
				//}
			}
		};

		if (det % 2 == 0) {
			//std::cout << "evennnnnnnnnnnnnnn" << std::endl;

			is_vec pivot;
			for (size_t i = 0; i < n; ++i) {
				pivot = adjA.col(i);
				int_type adj_col_2pow = vec_max2pow(pivot);
				if ((static_cast<int_type>(det) / adj_col_2pow) % 2 == 0) {
					pivot /= adj_col_2pow;
					//out() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					break;
				}

			}


			pivot = pivot.unaryExpr(std::ref(mod2_unaryop)).eval();

			/*for (size_t i = 0; i < pivot.size(); ++i)
				if (pivot[i] != 0 && pivot[i] != 1) {
					cout << "pivot[i]" << pivot[i] << endl;
					cin.ignore();
				}*/

			for (typename is_vec::InnerIterator it(pivot); it; ++it) {
				if (it.value() > 0)
					process_new_candidate(it.index(), pivot, 2);
			}
			/*for (size_t i = 0; i < n; ++i) {
				if (pivot[i] != 0)
					process_new_candidate(i, pivot, 2);
			}*/

		}
		else {
			//std::cout << "oddddddddddddddddodddd" << std::endl;

			is_vec pivot;
			//out() << A * adjA << std::endl;
			for (size_t i = 0; i < n; ++i) {
				pivot = adjA.col(i).unaryExpr(mod_unaryop(det));
				//out() << pivot.transpose() << std::endl;
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

			//out() << pivot.transpose() << std::endl;

			//std::cout << pivot << std::endl;
			for (typename is_vec::InnerIterator it(pivot); it; ++it) {
				if (it.value() > 0) {
					int_type lambda = _mod_equation_solve(it.value(), det - gcd(det, it.value()), det);
					is_vec t = (pivot * lambda).unaryExpr(mod_unaryop(det));

					//out() << t.transpose() << std::endl;

					process_new_candidate(it.index(), t, det);
				}
			}

			/*for (size_t i = 0; i < n; ++i) {
				if (pivot[i] > 0) {
					int_type lambda = _mod_equation_solve(pivot[i], det - gcd(det, pivot[i]), det);
					int_vec t = (pivot.array() * lambda).unaryExpr(mod_unaryop(det));

					process_new_candidate(i, t, det);
				}
			}*/
		}

		return std::move(result);
	}

	template<typename intT, typename ratT>
	typename ConicCutsGenerator<intT, ratT>::int_type 
		ConicCutsGenerator<intT, ratT>::_mod_equation_solve(int_type a, int_type b, int_type mod) const {
#ifdef _DEBUG
		assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif
		
		int_type x = 0, y = 0;
		int_type gcd_a_mod = 0; // = Internal::euclid_bezout(a, mod, x, y);

		mpz_gcdext(gcd_a_mod.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

		if (x < 0) x += mod;

#ifdef _DEBUG
		assert(x > 0);
		assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
		x *= (b / gcd_a_mod);

		return x;
	}

	template<typename intT, typename ratT>
	typename ConicCutsGenerator<intT, ratT>::int_type 
		ConicCutsGenerator<intT, ratT>::vec_max2pow(const is_vec& vec) const {

		int_type res = 0;

		mp_bitcnt_t min_bit_index = std::numeric_limits<mp_bitcnt_t>::max();
		for (typename is_vec::InnerIterator it(vec); it; ++it) {
			if (it.value() != 0) {
				mp_bitcnt_t bit_index = mpz_scan1(it.value().get_mpz_t(), 0);
				
					//std::cout << it.index() <<", " << it.value() << std::endl;
				min_bit_index = std::min(min_bit_index, bit_index);
			}
		}
		/*while (
			(vec_copy.unaryExpr(std::ref(mod2_unaryop)) == 0).all()
			) {
			res *= 2;
			vec_copy /= 2;
		}*/

		mpz_setbit(res.get_mpz_t(), min_bit_index);

		return std::move(res);
	}
}