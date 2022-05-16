#pragma once

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

#include <vector>
#include <numeric>
#include <functional>
#include <sstream>
#include <random>
#include <chrono>
#include <atomic>
#include <memory>
#include <fstream>

#include "gmp_eigen_supp.h"

namespace ILPCuts
{

	using namespace Eigen;

	template <typename T>
	class Atomic
	{
		std::atomic<T> _a;

	public:
		Atomic()
			: _a()
		{
		}

		Atomic(const std::atomic<T> &a)
			: _a(a.load())
		{
		}

		Atomic(const Atomic &other)
			: _a(other._a.load())
		{
		}

		Atomic &operator=(const Atomic &other)
		{
			_a.store(other._a.load());
		}

		Atomic &operator++()
		{
			++_a;
			return *this;
		}

		// Define prefix decrement operator.
		Atomic &operator--()
		{
			--_a;
			return *this;
		}

		bool operator<(T &a) const
		{
			return _a < a;
		}

		bool operator>(T &a) const
		{
			return _a > a;
		}

		bool operator==(T &a) const
		{
			return _a == a;
		}

		bool operator!=(T &a) const
		{
			return _a != a;
		}

		bool operator>=(T &a) const
		{
			return _a >= a;
		}
		bool operator<=(T &a) const
		{
			return _a <= a;
		}
	};

	class ToleranceErrorException : public std::runtime_error
	{
	public:
		ToleranceErrorException()
			: std::runtime_error("") {}

		ToleranceErrorException(std::string msg)
			: std::runtime_error(msg) {}
	};

	namespace Utility
	{

		template <typename inT, typename outT>
		struct ToSingle
		{
			outT operator()(const inT &input) const
			{
				return outT(input);
			}
		};

		template <>
		struct ToSingle<mpz_class, int64_t>
		{
			int64_t operator()(const mpz_class &input) const
			{
				return input.get_si();
			}
		};

		template <>
		struct ToSingle<mpz_class, int32_t>
		{
			int32_t operator()(const mpz_class &input) const
			{
				return input.get_si();
			}
		};

		template <>
		struct ToSingle<mpf_class, double>
		{
			double operator()(const mpf_class &input) const
			{
				return input.get_si();
			}
		};

		template <>
		struct ToSingle<mpq_class, double>
		{
			double operator()(const mpq_class &input) const
			{
				return input.get_d();
			}
		};

		template <typename intT, typename ratT>
		intT rat_to_int(const ratT &rat_x, const ratT &err_tolerance)
		{
			intT int_x = static_cast<intT>(round(rat_x));
			if (abs(rat_x - int_x) <= err_tolerance)
			{
				return int_x;
			}
			else
			{
				std::stringstream msg;
				msg << "rat_to_int tolerance error (" << abs(rat_x - int_x) << ")";
				throw ToleranceErrorException(msg.str());
			}
			return intT(0);
		}

		template <typename intT>
		intT mod2_unaryop(const intT &x)
		{
			return abs(x % intT(2));
		}

		template <typename intT>
		class mod_unaryop
		{
		private:
			const intT &_mod;

		public:
			mod_unaryop(const intT &mod)
				: _mod(mod) {}

			intT operator()(const intT &x) const
			{
				intT res(x % _mod);
				if (res < 0)
					res += _mod;
				return res;
			}
		};

		template <typename ratT>
		class Round
		{
		public:
			ratT operator()(const ratT &rat_x) const
			{
				return std::move(round(rat_x));
			}
		};

		template <typename ratT>
		class Floor
		{
		public:
			ratT operator()(const ratT &rat_x) const
			{
				return std::move(floor(rat_x));
			}
		};

		template <typename intT, typename ratT>
		class RatToInt_Unary
		{
		private:
			ratT _err_tolerance;

		public:
			RatToInt_Unary(const ratT &err_tolerance)
				: _err_tolerance(err_tolerance) {}

			intT operator()(const ratT &rat_x) const
			{
				return rat_to_int<intT, ratT>(rat_x, _err_tolerance);
			}
		};

		template <typename vecT>
		class SparseVecComp
		{
		public:
			bool operator()(const vecT &v1, const vecT &v2)
			{
				typename vecT::InnerIterator it1(v1);
				typename vecT::InnerIterator it2(v2);
				while (it1 && it2)
				{
					if (it1.index() > it2.index())
					{
						if (it1.value() < 0)
							return true;
						else
							return false;
					}
					else if (it1.index() < it2.index())
					{
						if (it2.value() > 0)
							return true;
						else
							return false;
					}
					else
					{
						if (it1.value() < it2.value())
							return true;
						else if (it1.value() > it2.value())
							return false;
					}

					++it1;
					++it2;
				}

				if (it1)
				{
					if (it1.value() < 0)
						return true;
					else
						return false;
				}

				if (it2)
				{
					if (it2.value() > 0)
						return true;
					else
						return false;
				}

				return false;
			}
		};

		template <typename vecT>
		class SparseVecEquiv
		{
		public:
			bool operator()(const vecT &v1, const vecT &v2)
			{
				typename vecT::InnerIterator it1(v1);
				typename vecT::InnerIterator it2(v2);

				while (it1 && it2)
				{

					if (it1.index() != it2.index())
						return false;
					else
					{
						if (it1.value() != it2.value())
							return false;
					}

					++it1;
					++it2;
				}

				return true;
			}
		};

		template <typename ratMatrixT, typename solverT>
		ratMatrixT compute_inverse(const ratMatrixT &A, solverT &solver);

		template <typename rationalT, typename solverT>
		SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT> &A, solverT &solver);

		template <typename rationalT>
		std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A);

		template <typename rationalT>
		rationalT compute_abs_det(const SparseMatrix<rationalT> &A);

		template <typename T>
		T euclid_bezout(
			const T &a,
			const T &b,
			T &u,
			T &v)
		{
			u = 1;
			v = 0;
			T u_new = v;
			T v_new = u;
			T r = abs(a);
			T r_new = abs(b);

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

			return r;
		}

		template <typename intT>
		intT mod_equation_solve(const intT &a, const intT &b, const intT &mod);

		template <>
		mpz_class mod_equation_solve(const mpz_class &a, const mpz_class &b, const mpz_class &mod)
		{
#ifdef _DEBUG
			assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif

			mpz_class x = 0, y = 0;
			mpz_class gcd_a_mod = 0; // = Internal::euclid_bezout(a, mod, x, y);

			mpz_gcdext(gcd_a_mod.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

			if (x < 0)
				x += mod;

#ifdef _DEBUG
			assert(x > 0);
			assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
			x *= (b / gcd_a_mod);

			return x;
		}

		template <typename rationalT>
		SparseVector<rationalT> system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b);

		template <typename rationalT, typename solverT>
		SparseVector<rationalT> system_solve(const SparseVector<rationalT> &b, solverT &solver);

		template <typename intT, typename randGenT>
		SparseVector<intT> bernouly_vec_generate(size_t n, randGenT &rgen);

		template <typename intT>
		class KeepIntFunc
		{
		private:
			const intT zero;

		public:
			KeepIntFunc() : zero(intT(0)) {}

			bool operator()(const intT &x) const { return x != zero; }
		};

		template <typename intT>
		void prune_zeroes_int(SparseMatrix<intT> &A)
		{
			A.prune(0, 1 / 2);
		}

		template <typename intT>
		void prune_zeroes_int(SparseVector<intT> &A)
		{
			A.prune(0, 1 / 2);
		}

		template <typename valT>
		valT max_norm(const SparseMatrix<valT> &A)
		{
			using value_type = valT;
			using mat_type = SparseMatrix<valT>;

			valT res(0);

			for (size_t k = 0; k < A.outerSize(); ++k)
				for (typename mat_type::InnerIterator it(A, k); it; ++it)
				{
					res = max(res, static_cast<valT>(abs(it.value())));
				}

			return std::move(res);
		}

		template <typename intT>
		intT vec_max2pow(const SparseVector<intT> &vec);

		template <>
		mpz_class vec_max2pow(const SparseVector<mpz_class> &vec)
		{
			using int_type = mpz_class;
			using is_vec = SparseVector<mpz_class>;

			int_type res = 0;

			mp_bitcnt_t min_bit_index = std::numeric_limits<mp_bitcnt_t>::max();
			for (typename is_vec::InnerIterator it(vec); it; ++it)
			{
				if (it.value() != 0)
				{
					mp_bitcnt_t bit_index = mpz_scan1(it.value().get_mpz_t(), 0);

					// std::cout << it.index() <<", " << it.value() << std::endl;
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

		template <typename intT>
		struct sparce_transform
		{
			using int_type = intT;

			size_t i;
			SparseVector<int_type> t;
			int_type denum;
		};

		template <typename int_type, typename rational_type>
		SparseVector<rational_type> compute_transformed_coeffs(const SparseVector<rational_type> &x,
															   const sparce_transform<int_type> &trans, double tolerance)
		{

			using rs_vec = SparseVector<rational_type>;
			using is_vec = SparseVector<int_type>;

			size_t n = x.rows();

			/*	cout << "N:   " << n << endl;
				cin.ignore();*/

			rs_vec res(n);

			int_type ti = trans.t.coeff(trans.i);
			rational_type xi = x.coeff(trans.i);

			if (xi == 0)
			{
				rs_vec xx = x * ti;
				res = xx / trans.denum;
			}
			else
			{
				typename rs_vec::InnerIterator x_it(x);
				typename is_vec::InnerIterator t_it(trans.t);
				bool negative_solution_found = false;

				while (t_it && x_it && !negative_solution_found)
				{
					if (t_it.index() < x_it.index())
					{
						negative_solution_found = true;
					}
					else if (x_it.index() < t_it.index())
					{
						rational_type y(x_it.value() * ti / trans.denum);
						res.insert(x_it.index()) = std::move(y);
						++x_it;
					}
					else
					{
						rational_type y((x_it.value() * ti - static_cast<rational_type>(t_it.value()) * xi) / trans.denum);
						if (y < rational_type(-tolerance))
							negative_solution_found = true;
						else
							res.insert(t_it.index()) = std::move(y);
						/*if (y > tolerance)
							res.insert(t_it.index()) = std::move(y);
						else if (y < -tolerance)
							negative_solution_found = true;*/
						++x_it;
						++t_it;
					}
				}

				if (!negative_solution_found)
				{
					if (t_it)
					{
						negative_solution_found = true;
					}
					else
					{
						while (x_it)
						{
							rational_type y(x_it.value() * ti / trans.denum);
							res.insert(x_it.index()) = std::move(y);
							++x_it;
						}

						res.coeffRef(trans.i) = xi;
					}
				}

				if (negative_solution_found)
					res.resize(0);
			}

			return std::move(res);
		}

		template <typename int_type>
		SparseMatrix<int_type> compute_adj_Q(size_t n, const sparce_transform<int_type> &trans)
		{
			using is_mat = SparseMatrix<int_type>;
			using is_vec = SparseVector<int_type>;

			const auto &t = trans.t;
			size_t i = trans.i;
			const auto &denum = trans.denum;
			const auto &ti = t.coeff(i);

			is_mat MM(n, n);
			using Tripl = Triplet<int_type>;
			std::vector<Tripl> triplets;
			for (size_t j = 0; j < n; ++j)
			{
				if (j != i)
					triplets.push_back(Tripl(j, j, ti));
				else
					triplets.push_back(Tripl(j, j, denum));
			}

			for (typename is_vec::InnerIterator it(t); it; ++it)
			{
				if (it.value() > 0 && it.index() != i)
					triplets.emplace_back(Tripl(
						it.index(),
						i,
						int_type(-it.value())));
			}

			MM.setFromTriplets(triplets.begin(), triplets.end());

			return std::move(MM);
		}

		// class PrimesDownloader
		// {
		// private:
		// 	std::ifstream file;

		// public:
		// 	static std::string get_default_primes_base_path()
		// 	{
		// 		return "../primes_table/primes_mill.txt";
		// 	}

		// 	PrimesDownloader(std::string path = get_default_primes_base_path())
		// 		: file(path) {}

		// };

		namespace primes
		{
			const std::string default_path = "../primes_table/primes_mill.txt";
			const size_t max_download = 1000000;

			template <typename ptype>
			std::vector<ptype> download(size_t max_num, std::string path = default_path)
			{
				std::ifstream file(path);
				std::vector<ptype> res;

				size_t i = 0;
				while (file && (i++ < max_num))
				{
					res.push_back(0);
					file >> res.back();
				}

				return std::move(res);
			}

			template <typename int_type, typename lint_type>
			int_type find_divisor(const lint_type &val, const std::vector<int_type> &primes_data)
			{
				for (int_type p : primes_data)
				{
					if (val % p == 0)
						return p;
				}

				return 0;
			}
		}

	};

	template <typename rationalT, typename solverT>
	SparseMatrix<rationalT> Utility::compute_inverse(const SparseMatrix<rationalT> &A, solverT &solver)
	{
		using rat_mat = SparseMatrix<rationalT>;

		size_t n = A.cols();

		rat_mat I(n, n);
		I.reserve(VectorXi::Constant(n, 2));
		for (size_t i = 0; i < n; ++i)
			I.insert(i, i) = 1;

		rat_mat res(n, n);
		/*for (size_t i = 0; i < n; ++i) {
			res.col(i) = solver.solve(I.col(i));
		}*/
		res = std::move(solver.solve(I));

		return std::move(res);
	}

	template <typename rationalT>
	std::tuple<SparseMatrix<rationalT>, rationalT> Utility::compute_abs_adjugate_and_det(const SparseMatrix<rationalT> &A)
	{

		using rat_type = rationalT;
		using rat_mat = SparseMatrix<rat_type>;

		SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

		rat_type det = abs(decomp_LU.determinant());
		rat_mat adjA = compute_inverse<rat_type>(A, decomp_LU) * det;

		return std::make_tuple(std::move(adjA), std::move(det));
	}

	template <typename rationalT>
	rationalT Utility::compute_abs_det(const SparseMatrix<rationalT> &A)
	{
		using rat_type = rationalT;
		using rat_mat = SparseMatrix<rat_type>;

		SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

		rat_type det = abs(decomp_LU.determinant());

		return std::move(det);
	}

	template <typename rationalT>
	SparseVector<rationalT> Utility::system_solve(const SparseMatrix<rationalT> &A, const SparseVector<rationalT> &b)
	{
		using rat_type = rationalT;
		using rat_mat = SparseMatrix<rat_type>;
		using rat_vec = SparseVector<rationalT>;

		SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

		return std::move(system_solve(b, decomp_LU));
	}

	template <typename rationalT, typename solverT>
	SparseVector<rationalT> Utility::system_solve(const SparseVector<rationalT> &b, solverT &solver)
	{
		using rat_type = rationalT;
		using rat_vec = SparseVector<rationalT>;
		using solver_type = solverT;

		return std::move(solver.solve(b));
	}

	template <typename intT, typename randGenT>
	SparseVector<intT> Utility::bernouly_vec_generate(size_t n, randGenT &rgen)
	{
		using int_type = intT;
		using is_vec = SparseVector<int_type>;

		std::bernoulli_distribution bd(0.5);

		// generation of a random rhs
		is_vec res(n);
		res.reserve((2 * n) / 3);
		for (size_t i = 0; i < n; ++i)
		{
			bool flip = bd(rgen);
			if (flip)
				res.insert(i) = 1.0;
		}

		return std::move(res);
	}

}
