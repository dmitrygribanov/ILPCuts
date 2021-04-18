#pragma once

//#include<gcd.hpp>

#include<Eigen\SparseCore>
#include<Eigen\Dense>

#include<vector>
#include<numeric>
#include<functional>
#include<sstream>

namespace ILPCuts {

	using namespace Eigen;

	namespace Internal {
		template <typename T>
		T euclid_bezout
		(
			const T& a,
			const T& b,
			T& u,
			T& v
		)
		{
			u = 1;
			v = 0;
			T u_new = v;
			T v_new = u;
			T r = std::abs(a);
			T r_new = std::abs(b);

			while (r_new != 0)
			{
				T q = r / r_new;

				T aux = u_new;
				u_new = u - q * u_new;
				u = aux;

				aux = v_new;
				v_new = v - q * v_new;
				v = aux;

				aux = r_new;
				r_new = r - q * r_new;
				r = aux;
			}

			/*if (is_negative(a))
				opposite(&u);
			if (is_negative(b))
				opposite(&v);*/

				//_Internal::euclid_norm(r, u, v);

				//ARAGELI_ASSERT_1(a * u + b * v == r);

			return r;
		}

		/*template <typename T>
	T euclid(T a, T b)
	{
		while ((a %= b) != 0)
		{
			std::swap(a, b);
		}
		return b;
	}

	template <typename T>
	inline T gcd
	(
		const T& a,
		const T& b
	)
	{
		return euclid(a, b);
	}*/
	}

	template<typename intT, typename rationalT,
		template<typename> class matrixT,
		template<typename> class vectorT>
	class ConicCutsGenerator_Eigen {
	public:

		using rational_type = rationalT;
		using int_type = intT;

		using int_mat = matrixT<int_type>;
		using rat_mat = matrixT<rational_type>;

		using int_vec = vectorT<int_type>;
		using rat_vec = vectorT<rational_type>;

	private:

		std::stringstream _out;

		int_type vec_gcd(const int_vec& vec) const {
			int_type res = std::gcd(vec[0], vec[1]);

			for (size_t i = 2; i < vec.size(); ++i) {
				res = std::gcd(res, vec[i]);
			}

			return res;
		}

		void normalize(int_vec & v) const {
			int_type gcd = vec_gcd(v);
			v = v.array() / gcd;
		}

		int_type vec_max2pow(const int_vec& vec) const {
			int_type res = 1;

			int_vec vec_copy = vec;

			while ( 
				(vec_copy.array().unaryExpr(std::ref(mod2_unaryop)) == 0).all()
				) {
				res *= 2;
				vec_copy /= 2;
			}

			return res;
		}

		static int_type mod2_unaryop(const int_type& x) {
			return x % 2;
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

		int_type _mod_equation_solve(int_type a, int_type b, int_type mod) const {
#ifdef _DEBUG
			assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif
			int_type x = 0, y = 0;
			int_type gcd_a_mod = Internal::euclid_bezout(a, mod, x, y);

			if (x < 0) x += mod;

#ifdef _DEBUG
			assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
			x *= (b / gcd_a_mod);

			return x;
		}

		static int_type rat_to_int(const rational_type& rat_x, const rational_type& err_tolerance) {
			int_type int_x = static_cast<int_type>(std::round(rat_x));
			if (std::abs(rat_x - int_x) <= err_tolerance) {
				return int_x;
			}
			else {
				std::stringstream msg;
				msg << "rat_to_int tolerance error (" << std::abs(rat_x - int_x) << ")";
				throw ToleranceErrorException(msg.str());
			}
		}

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

		struct cone_traits {
			int_mat A;
			int_mat adjA;
			int_type det;
		};

		std::vector<cone_traits> compute_layer(const cone_traits& cone, const rat_mat& directions_pool,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {

			std::vector<cone_traits> result;
			const int_mat& A = cone.A;
			const int_mat& adjA = cone.adjA;
			const int_type& det = cone.det;
			size_t n = A.cols();

			/*void process_new_candidate(const int_mat & a, const rat_vec & t) {
				result.size();
			};*/

			auto process_new_candidate = [&](size_t i, const int_vec& t, const int_type& div)
			{
				cone_traits new_cone;

				try {
					new_cone.A = A;
					new_cone.A.col(i) = ((A * t).cast<rational_type>().array() / div)
						.unaryExpr(RatToInt_Unary(err_tolerance));
					//.round()
					//.cast<int_type>();
				}
				catch (ToleranceErrorException exc) {
					out() << "b-column construction tolerance error: " << exc.what() << '\n';
				}

				PartialPivLU<rat_mat> decomp_LU(new_cone.A.cast<rational_type>());

				try {

					new_cone.det = rat_to_int(std::abs(decomp_LU.determinant()), err_tolerance);

				}
				catch (ToleranceErrorException exc) {
					out() << "Determinant computing tolerance error: " << exc.what() << '\n';
				}

				try {
					new_cone.adjA = (decomp_LU.inverse().array() * new_cone.det)
						.round()
						.cast<int_type>();
				}
				catch (ToleranceErrorException exc) {
					out() << "Inversion tolerance error: " << exc.what() << '\n';
				}

				for (size_t i = 0; i < directions_pool.cols(); ++i) {
					rat_vec duals = new_cone.adjA.cast<rational_type>() * directions_pool.col(i);
					if (((duals.array() + tolerance) > 0).all()) {
						result.push_back(std::move(new_cone));
					}
				}
			};

			if (det % 2 == 0) {
				int_vec pivot;
				for (size_t i = 0; i < n; ++i) {
					pivot = adjA.col(i);
					int_type adj_col_2pow = vec_max2pow(pivot);
					if ((static_cast<int_type>(det) / adj_col_2pow) % 2 == 0) {
						pivot /= adj_col_2pow;
						break;
					}
				}

				pivot = pivot.unaryExpr(std::ref(mod2_unaryop)).eval();

				for (size_t i = 0; i < n; ++i) {
					process_new_candidate(i, pivot, 2);
				}

			}
			else {
				int_vec pivot = adjA.array().unaryExpr(mod_unaryop(det));
				for (size_t i = 0; i < n; ++i) {
					if (pivot[i] > 0) {
						int_type lambda = _mod_equation_solve(pivot[i], det - std::gcd(det, pivot[i]), det);
						int_vec t = (pivot.array() * lambda).unaryExpr(mod_unaryop(det));

						process_new_candidate(i, t, det);
					}
				}
			}

			return std::move(result);
		}

		std::vector<int_mat> partial_uni_decomp
		(const int_mat& A, const rat_mat& directions_pool,
			const int_mat& adjA, const int_type& det,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {

			std::vector<cone_traits> current_layer;
			std::vector<cone_traits> next_layer;
			std::vector<int_mat> result;

			current_layer.push_back(std::move(cone_traits{ A,adjA,det }));

			while (!current_layer.empty()) {

				for (const cone_traits& node : current_layer) {
					std::vector<cone_traits> next_node_layer = std::move(compute_layer(node, directions_pool, tolerance, err_tolerance));
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

	public:

		std::stringstream& out() { return _out; }

		class DirectionsStrategyBase {
		public:
			virtual rat_mat generate_directons() = 0; // { return rat_mat();  };
		};

		class GreadyDirection :public DirectionsStrategyBase {
		private:
			const rat_vec& _c;

		public:
			GreadyDirection(const rat_vec& c)
				:_c(c) {}

			rat_mat generate_directons() {
				return std::move(rat_mat(_c));
			}
		};

		std::vector<int_vec> generate_cuts(const int_mat& A, const int_mat& adjA, const rat_vec& optLP,
			DirectionsStrategyBase & directions_strategy,
			const int_type& det,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {
			
			std::vector<int_vec> res_cuts;

			rat_mat directions_pool = std::move( directions_strategy.generate_directons() );

			std::vector<int_mat> uni_mat_pool = std::move( partial_uni_decomp(A, directions_pool,
				adjA, det, tolerance, err_tolerance) );
			
			for (const int_mat& uni_mat : uni_mat_pool) {
				rat_vec rat_b = uni_mat.cast<rational_type>() * optLP;
				int_vec rounded_b = rat_b.array().round().cast<int_type>();
				int_vec floor_b = rat_b.array().floor().cast<int_type>();

				rat_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).array().abs();

				for (size_t i = 0; i < diff_b.size(); ++i) {
					if (diff_b[i] > tolerance) {
						int_vec cut = uni_mat.row(i);
						cut << floor_b[i];
						res_cuts.push_back(std::move( cut ));
					}
				}
			}

			return std::move(res_cuts);
		}

		std::vector<int_vec> generate_cuts(const int_mat& A, const rat_vec& optLP,
			DirectionsStrategyBase& directions_strategy,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {

			PartialPivLU<rat_mat> decomp_LU(A.cast<rational_type>());

			int_type det = 0;

			try {
				det = rat_to_int(std::abs(decomp_LU.determinant()), err_tolerance);
			}
			catch (ToleranceErrorException exc) {
				out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			}

			try {
				int_mat adjA = (decomp_LU.inverse().array() * det).unaryExpr(RatToInt_Unary(err_tolerance));

				return std::move(generate_cuts(A, adjA, optLP, directions_strategy,
					det, tolerance, err_tolerance));
			}
			catch (ToleranceErrorException exc) {
				out() << "Root: inverse tolerance error: " << exc.what() << '\n';
			}

		}

		class ToleranceErrorException : public std::exception {
		public:
			ToleranceErrorException()
				:std::exception() {}

			ToleranceErrorException(std::string msg)
				:std::exception(msg.c_str()) {}
		};
	};

	template<typename T>
	using DenceMat = Matrix<T, Dynamic, Dynamic>;

	template<typename T>
	using DenceVec = Matrix<T, Dynamic, 1>;

	template<typename T>
	using SparceMat = SparseMatrix<T>;

	template<typename T>
	using SparceVec = SparseVector<T>;

	template<typename intT = long long int, typename ratT = double>
	using ConicCutsGenerator_EigenDence = ConicCutsGenerator_Eigen<intT, ratT, 
		DenceMat, DenceVec >;

	template<typename intT = long long int, typename ratT = double>
	using ConicCutsGenerator_EigenSparce = ConicCutsGenerator_Eigen<intT, ratT,
		SparceMat, SparceVec >;
}
