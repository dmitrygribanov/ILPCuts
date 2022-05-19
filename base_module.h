#pragma once

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

#include <vector>
#include <numeric>
#include <functional>
#include <random>
#include <chrono>
#include <atomic>
#include <memory>

#include "gmp_eigen_supp.h"

namespace ILPCuts
{
	using namespace Eigen;
	using namespace std;
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
			const double &_err_tolerance;

		public:
			explicit RatToInt_Unary(const double &err_tolerance)
				: _err_tolerance(err_tolerance) {}

			intT operator()(const ratT &rat_x) const
			{
				return rat_to_int<intT, ratT>(rat_x, _err_tolerance);
			}

			const double &get_err_tolerance() const { return _err_tolerance; }
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

		template <>
		int64_t mod_equation_solve(const int64_t &a, const int64_t &b, const int64_t &mod)
		{
#ifdef _DEBUG
			assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif

			int64_t x = 0, y = 0;
			int64_t gcd_a_mod = euclid_bezout(a, mod, x, y);

			if (x < 0)
				x += mod;

#ifdef _DEBUG
			assert(x > 0);
			assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
			x *= (b / gcd_a_mod);

			return x;
		}

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
			A.prune(0, 0.5);
		}

		template <typename intT>
		void prune_zeroes_int(SparseVector<intT> &A)
		{
			A.prune(0, 0.5);
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

	};

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
				res.insert(i) = 1;
		}

		return std::move(res);
	}

}
