#pragma once

//#include<gcd.hpp>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>

#include <vector>
#include <numeric>
#include <functional>
#include <sstream>
#include <random>

#include "base_module.h"

namespace ILPCuts
{

	using namespace Eigen;

	template <typename intT, typename rationalT, template <typename> class matrixT, template <typename> class vectorT>
	class ConicCutsGenerator_Eigen
	{
	public:
		using rational_type = rationalT;
		using int_type = intT;

		using int_mat = matrixT<int_type>;
		using rat_mat = matrixT<rational_type>;

		using int_vec = vectorT<int_type>;
		using rat_vec = vectorT<rational_type>;

		/*typedef rationalT rational_type;
		typedef intT int_type;

		typedef matrixT<int_type> int_mat;
		typedef matrixT<rational_type> rat_mat;

		typedef vectorT<int_type> int_vec;
		typedef vectorT<rational_type> rat_vec;*/

	private:
		std::ostream *_out;
		std::ostream *_err;

		int_type vec_gcd(const int_vec &vec) const
		{
			int_type res = gcd(vec[0], vec[1]);

			for (size_t i = 2; i < vec.size(); ++i)
			{
				res = gcd(res, vec[i]);
			}

			return res;
		}

		void normalize(int_vec &v) const
		{
			int_type gcd = vec_gcd(v);
			v = v.array() / gcd;
		}

		int_type vec_max2pow(const int_vec &vec) const
		{
			int_type res = 1;

			int_vec vec_copy = vec;

			while (
				(vec_copy.array().unaryExpr(std::ref(mod2_unaryop)) == 0).all())
			{
				res *= 2;
				vec_copy /= 2;
			}

			return res;
		}

		static int_type mod2_unaryop(const int_type &x)
		{
			return abs(x % 2);
		}

		class mod_unaryop
		{
		private:
			const int_type &_mod;

		public:
			mod_unaryop(const int_type &mod)
				: _mod(mod) {}

			int_type operator()(const int_type &x) const
			{
				int_type res = x % _mod;
				if (res < 0)
					res += _mod;
				return res;
			}
		};

		int_type _mod_equation_solve(int_type a, int_type b, int_type mod) const
		{
#ifdef _DEBUG
			assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif
			int_type x = 0, y = 0;
			int_type gcd_a_mod = Utility::euclid_bezout(a, mod, x, y);

			if (x < 0)
				x += mod;

#ifdef _DEBUG
			assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
			x *= (b / gcd_a_mod);

			return x;
		}

		static int_type rat_to_int(const rational_type &rat_x, const rational_type &err_tolerance)
		{
			// std::cout << rat_x << std::endl;
			int_type int_x = static_cast<int_type>(round(rat_x));
			/*std::cout << int_x << std::endl;
			std::cout << abs(rat_x - int_x) << std::endl;*/
			if (abs(rat_x - int_x) <= err_tolerance)
			{
				return int_x;
			}
			else
			{
				std::stringstream msg;
				msg << "rat_to_int tolerance error (" << abs(rat_x - int_x) << ")";
				std::cout << abs(rat_x - int_x) << std::endl;
				throw ToleranceErrorException(msg.str());
			}
			return 0;
		}

		class RatToInt_Unary
		{
		private:
			rational_type _err_tolerance;

		public:
			RatToInt_Unary(const rational_type &err_tolerance)
				: _err_tolerance(err_tolerance) {}

			int_type operator()(const rational_type &rat_x) const
			{
				return rat_to_int(rat_x, _err_tolerance);
			}
		};

		struct cone_traits
		{
			int_mat A;
			int_mat adjA;
			int_type det;
		};

		std::vector<cone_traits> compute_layer(const cone_traits &cone, const rat_mat &directions_pool,
											   std::vector<size_t> &directions_evaluation_counter,
											   const rational_type &tolerance = 0.00001,
											   const rational_type &err_tolerance = 0.0001)
		{

			std::vector<cone_traits> result;
			const int_mat &A = cone.A;
			const int_mat &adjA = cone.adjA;
			const int_type &det = cone.det;
			size_t n = A.cols();

			std::cout << "current determinant: " << det << std::endl;

			auto process_new_candidate = [&](size_t i, const int_vec &t, const int_type &div)
			{
				cone_traits new_cone;

				try
				{
					new_cone.A = A;
					new_cone.A.col(i) = ((A * t).cast<rational_type>().array() / div)
											.unaryExpr(RatToInt_Unary(err_tolerance));

					// std::cout << A << std::endl;
					// std::cout << new_cone.A << std::endl;
					//.round()
					//.cast<int_type>();
				}
				catch (ToleranceErrorException exc)
				{
					out() << "b-column construction tolerance error: " << exc.what() << '\n';
					err() << "b-column construction tolerance error: " << exc.what() << '\n';
				}

				PartialPivLU<rat_mat> decomp_LU(new_cone.A.cast<rational_type>());

				try
				{

					new_cone.det = rat_to_int(abs(decomp_LU.determinant()), err_tolerance);
					/*if (det % 2 == 0) {
						cout << "old det:" << det << endl;
						cout << "new det:" << new_cone.det << endl;
					}*/
				}
				catch (ToleranceErrorException exc)
				{
					out() << "Determinant computing tolerance error: " << exc.what() << '\n';
					err() << "Determinant computing tolerance error: " << exc.what() << '\n';
				}

				try
				{
					new_cone.adjA = (decomp_LU.inverse().array() * new_cone.det)
										.round()
										.cast<int_type>();
				}
				catch (ToleranceErrorException exc)
				{
					out() << "Inversion tolerance error: " << exc.what() << '\n';
					err() << "Inversion tolerance error: " << exc.what() << '\n';
				}

				for (size_t i = 0; i < directions_pool.cols(); ++i)
				{
					// cout << "direction: " << directions_pool.col(i).transpose() << endl;
					rat_vec duals = new_cone.adjA.cast<rational_type>() * directions_pool.col(i);
					// rat_vec duals = adjA.cast<rational_type>() * directions_pool.col(i);
					// cout << "duals: " << duals << endl;

					/*int_mat I = A * adjA;
					cout << I << endl;
					cin.ignore();*/

					if (((duals.array() + tolerance) > 0).all())
					{
						result.push_back(std::move(new_cone));
						std::cout << "adding step done" << std::endl;
						// cin.ignore();
					}
				}
			};

			if (det % 2 == 0)
			{
				std::cout << "evennnnnnnnnnnnnnn" << std::endl;

				int_vec pivot;
				for (size_t i = 0; i < n; ++i)
				{
					pivot = adjA.col(i);
					int_type adj_col_2pow = vec_max2pow(pivot);
					if ((static_cast<int_type>(det) / adj_col_2pow) % 2 == 0)
					{
						pivot /= adj_col_2pow;
						break;
					}
				}

				pivot = pivot.unaryExpr(std::ref(mod2_unaryop)).eval();

				/*for (size_t i = 0; i < pivot.size(); ++i)
					if (pivot[i] != 0 && pivot[i] != 1) {
						cout << "pivot[i]" << pivot[i] << endl;
						cin.ignore();
					}*/

				for (size_t i = 0; i < n; ++i)
				{
					if (pivot[i] != 0)
						process_new_candidate(i, pivot, 2);
				}
			}
			else
			{
				std::cout << "oddddddddddddddddodddd" << std::endl;

				int_vec pivot;
				for (size_t i = 0; i < n; ++i)
				{
					pivot = adjA.col(i).unaryExpr(mod_unaryop(det));
					if ((pivot.array() != 0).any())
						break;
				}

				// std::cout << pivot << std::endl;
				for (size_t i = 0; i < n; ++i)
				{
					if (pivot[i] > 0)
					{
						int_type lambda = _mod_equation_solve(pivot[i], det - gcd(det, pivot[i]), det);
						int_vec t = (pivot.array() * lambda).unaryExpr(mod_unaryop(det));

						process_new_candidate(i, t, det);
					}
				}
			}

			return std::move(result);
		}

		// static int_type foo(int_type) { return 0; }

		std::vector<int_mat> partial_uni_decomp(const int_mat &A, const rat_mat &directions_pool,
												const int_mat &adjA, const int_type &det,
												const rational_type &tolerance = 0.00001,
												const rational_type &err_tolerance = 0.0001)
		{

			std::vector<cone_traits> current_layer;
			std::vector<cone_traits> next_layer;
			std::vector<int_mat> result;

			current_layer.push_back(std::move(cone_traits{A, adjA, det}));

			while (!current_layer.empty())
			{

				for (const cone_traits &node : current_layer)
				{
					std::vector<cone_traits> next_node_layer = std::move(compute_layer(node, directions_pool, tolerance, err_tolerance));
					for (cone_traits &next_node : next_node_layer)
					{
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
		ConicCutsGenerator_Eigen(std::ostream *outs = &std::cout, std::ostream *errs = &std::cout)
			: _out(outs), _err(errs) {}

		std::ostream &out() { return *_out; }

		std::ostream &err() { return *_err; }

		class DirectionsStrategyBase
		{
		public:
			virtual rat_mat generate_directons() = 0; // { return rat_mat();  };
		};

		class GreadyDirection : public DirectionsStrategyBase
		{
		private:
			const rat_vec &_c;

		public:
			GreadyDirection(const rat_vec &c)
				: _c(c) {}

			rat_mat generate_directons()
			{
				return std::move(rat_mat(_c));
			}
		};

		class SomeDirection : public DirectionsStrategyBase
		{
		private:
			rat_vec _direction;

		public:
			SomeDirection(const int_mat &A)
			{
				size_t n = A.cols();
				_direction = A.cast<rational_type>().rowwise().sum().array() / n;
			}

			rat_mat generate_directons()
			{
				return std::move(rat_mat(_direction));
			}
		};

		template <typename randGenT>
		class ConvUniformAndGready : public GreadyDirection
		{
		public:
			using gen_type = randGenT;

		private:
			const int_mat &_A;
			const size_t _rand_dir_num;
			gen_type &_gen;

		public:
			ConvUniformAndGready(gen_type &gen, const rat_vec &c, const int_mat &A, size_t rand_dir_num = 3)
				: GreadyDirection(c), _gen(gen), _A(A), _rand_dir_num(rand_dir_num) {}

			rat_mat generate_directons()
			{
				rat_mat result_directions = std::move(GreadyDirection::generate_directons());
				size_t n = _A.cols();

				result_directions.conservative_resize(n, _rand_dir_num + 1);

				// rat_mat NA = _A.cast<rational_type>().colwise().normalize();
				std::uniform_real_distribution<rational_type> udistr(0.0, 1.0);

				for (size_t i = 0; i < _rand_dir_num; ++i)
				{
					rat_vec new_dir = _A.col(0).cast<rational_type>();
					for (size_t j = 1; j < n; ++j)
					{
						rational_type alpha = udistr(_gen);
						new_dir = (new_dir * alpha + _A.col(j).cast<rational_type>() * (1.0 - alpha)).eval();
					}
					result_directions.col(i + 1) = std::move(new_dir);
				}

				return std::move(result_directions);
			}
		};

		std::vector<int_vec> generate_cuts(const int_mat &A, const int_mat &adjA, const rat_vec &optLP,
										   DirectionsStrategyBase &directions_strategy,
										   const int_type &det,
										   const rational_type &tolerance = 0.00001,
										   const rational_type &err_tolerance = 0.0001)
		{

			std::vector<int_vec> res_cuts;

			rat_mat directions_pool = std::move(directions_strategy.generate_directons());

			std::vector<int_mat> uni_mat_pool = std::move(partial_uni_decomp(A, directions_pool,
																			 adjA, det, tolerance, err_tolerance));

			for (const int_mat &uni_mat : uni_mat_pool)
			{
				rat_vec rat_b = uni_mat.cast<rational_type>() * optLP;
				int_vec rounded_b = rat_b.array().round().cast<int_type>();
				int_vec floor_b = rat_b.array().floor().cast<int_type>();

				rat_vec diff_b = (rat_b - rounded_b.cast<rational_type>()).array().abs();

				/*std::cout << uni_mat << std::endl;
				std::cout << optLP.transpose() << std::endl;
				std::cout << rat_b.transpose() << std::endl;
				std::cout << rounded_b.transpose() << std::endl;
				std::cout << floor_b.transpose() << std::endl;
				std::cout << diff_b.transpose() << std::endl;*/

				for (size_t i = 0; i < diff_b.size(); ++i)
				{
					if (diff_b[i] > tolerance)
					{
						int_vec cut = uni_mat.row(i);
						cut.conservativeResize(cut.size() + 1);
						cut[cut.size() - 1] = floor_b[i];
						out() << "cut has been added: (" << cut.transpose() << ")\n";
						res_cuts.push_back(std::move(cut));
					}
				}
			}

			return std::move(res_cuts);
		}

		std::vector<int_vec> generate_cuts(const int_mat &A, const rat_vec &optLP,
										   DirectionsStrategyBase &directions_strategy,
										   const rational_type &tolerance = 0.00001,
										   const rational_type &err_tolerance = 0.0001)
		{

			PartialPivLU<rat_mat> decomp_LU(A.cast<rational_type>());

			int_type det = 0;

			// std::cout << A << std::endl;

			try
			{
				det = rat_to_int(abs(decomp_LU.determinant()), err_tolerance);
			}
			catch (ToleranceErrorException exc)
			{
				out() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
				err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			}

			try
			{
				int_mat adjA = (decomp_LU.inverse().array() * det).unaryExpr(RatToInt_Unary(err_tolerance));

				return std::move(generate_cuts(A, adjA, optLP, directions_strategy,
											   det, tolerance, err_tolerance));
			}
			catch (ToleranceErrorException exc)
			{
				out() << "Root: inverse tolerance error: " << exc.what() << '\n';
				err() << "Root: determinant computing tolerance error: " << exc.what() << '\n';
			}
		}

		class ToleranceErrorException : public std::exception
		{
		public:
			ToleranceErrorException()
				: std::exception() {}

			ToleranceErrorException(std::string msg)
				: std::exception(msg.c_str()) {}
		};
	};

	template <typename T>
	using DenceMat = Matrix<T, Dynamic, Dynamic>;

	template <typename T>
	using DenceVec = Matrix<T, Dynamic, 1>;

	template <typename T>
	using SparceMat = SparseMatrix<T>;

	template <typename T>
	using SparceVec = SparseVector<T>;

	template <typename intT = long long int, typename ratT = double>
	using ConicCutsGenerator_EigenDence = ConicCutsGenerator_Eigen<intT, ratT,
																   DenceMat, DenceVec>;

	template <typename intT = long long int, typename ratT = double>
	using ConicCutsGenerator_EigenSparce = ConicCutsGenerator_Eigen<intT, ratT,
																	SparceMat, SparceVec>;
}
