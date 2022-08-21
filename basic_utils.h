#pragma once

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

#include <vector>
#include <random>
#include <atomic>
#include <string>
#include <sstream>
#include <future>
#include <thread>
#include <functional>

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

	template <typename callableT>
	bool run_async_function_with_time_limit(double time_limit, callableT fun)
	{
		// create wrapper for a function to create a tread
		packaged_task<void()> task(fun);

		// taking shared information for the task called "future"
		auto task_future = task.get_future();

		// creating a tread
		thread task_thread(std::move(task));

		// waiting for the timeout
		auto wait_res = task_future.wait_for(chrono::duration<double>(time_limit));

		// if the work is not done, detach the thread
		// and return failure
		if (wait_res == future_status::timeout)
		{
			task_thread.detach();

			return false;
		}

		// if the work is done, call "future::get()" to rethrow exceptions
		// stored in shared state
		task_future.get();

		// waiting for a correct termination of the thread
		task_thread.join();

		return true;
	}

	class Exception : public std::runtime_error
	{
	public:
		Exception()
			: std::runtime_error("ILPCuts lib error") {}

		Exception(std::string msg)
			: std::runtime_error(msg) {}
	};

	class ToleranceErrorException : public Exception
	{
	public:
		ToleranceErrorException()
			: Exception("ILPCuts lib error: a tolerance-error has occurred.") {}

		ToleranceErrorException(std::string msg)
			: Exception(msg) {}
	};

	template <typename inT, typename outT>
	struct conv
	{
		static outT make(const inT &a)
		{
			return outT(a);
		}
	};

	template <typename inT, typename outVecT>
	struct conv<SparseVector<inT>, outVecT>
	{
		static outVecT make(const SparseVector<inT> &a)
		{
			using out_val_type = typename outVecT::value_type;

			size_t n = a.size();

			vector<out_val_type> res(n, out_val_type(0));

			for (typename SparseVector<inT>::InnerIterator it(a); it; ++it)
			{
				res[it.index()] = conv<inT, out_val_type>::make(it.value());
			}

			return res;
		}
	};

	template <typename inVecT, typename outT>
	struct conv<inVecT, SparseVector<outT>>
	{
		static SparseVector<outT> make(const inVecT &a)
		{
			using in_val_type = typename inVecT::value_type;

			size_t n = a.size();

			SparseVector<outT> res(n);
			for (size_t i = 0; i < n; ++i)
			{
				if (a[i] != 0)
					res.insertBack(i) = conv<in_val_type, outT>::make(a[i]);
			}

			return res;
		}
	};
	namespace Utility
	{

		template <typename intT, typename ratT>
		intT rat_to_int(const ratT &rat_x, double err_tolerance)
		{
			ratT rounded_x = round(rat_x);
			ratT diff = abs(rat_x - rounded_x);
			if (diff <= err_tolerance)
			{
				return conv<ratT, intT>::make(rounded_x);
			}
			else
			{
				std::stringstream msg;
				msg << "ILPCuts tolerance error:\n diff = (" << diff
					<< ")\n value = (" << rat_x
					<< ")\n rounded value = (" << rounded_x << ")";
				throw ToleranceErrorException(msg.str());
			}
		}

		template <typename intT>
		intT mod2_unaryop(const intT &x)
		{
			return abs(x % intT(2));
		}

		template <>
		int64_t mod2_unaryop(const int64_t &x)
		{
			return x & 1;
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
		class Round_Unary
		{
		public:
			ratT operator()(const ratT &rat_x) const
			{
				return round(rat_x);
			}
		};

		template <typename ratT>
		class Floor_Unary
		{
		public:
			ratT operator()(const ratT &rat_x) const
			{
				return floor(rat_x);
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

		template <typename T>
		class SparseVecLess
		{
		public:
			bool operator()(const SparseVector<T> &v1, const SparseVector<T> &v2) const
			{
				typename SparseVector<T>::InnerIterator it1(v1);
				typename SparseVector<T>::InnerIterator it2(v2);

				while (it1 && it2)
				{
					if (it1.index() > it2.index())
					{
						if (it2.value() > 0)
							return true;
						else
							return false;
					}
					else if (it1.index() < it2.index())
					{
						if (it1.value() < 0)
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

		template <typename T>
		class SparseVecEquiv
		{
		public:
			bool operator()(const SparseVector<T> &v1, const SparseVector<T> &v2) const
			{
				typename SparseVector<T>::InnerIterator it1(v1);
				typename SparseVector<T>::InnerIterator it2(v2);

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
		SparseVector<intT> uniform_vec_generate(size_t n, randGenT &rgen, const intT &a, const intT &b)
		{
			using int_type = intT;
			using is_vec = SparseVector<int_type>;

			std::uniform_int_distribution<int_type> ud(a, b - 1);

			// generation of a random rhs
			is_vec res(n);
			res.reserve(size_t(0.7 * n));
			for (size_t i = 0; i < n; ++i)
			{
				int_type val = ud(rgen);
				if (val != 0)
					res.insert(i) = val;
			}

			return res;
		}

		template <typename intT, typename randGenT>
		SparseVector<intT> bernouly_vec_generate(size_t n, randGenT &rgen)
		{
			using int_type = intT;
			using is_vec = SparseVector<int_type>;

			std::bernoulli_distribution bd(0.5);

			// generation of a random rhs
			is_vec res(n);
			res.reserve(size_t(0.7 * n));
			for (size_t i = 0; i < n; ++i)
			{
				bool flip = bd(rgen);
				if (flip)
					res.insert(i) = 1;
			}

			return res;
		}

		template <typename intT>
		class KeepIntFunc
		{
		private:
			const intT zero;

		public:
			KeepIntFunc() : zero(0) {}

			bool operator()(size_t i, size_t j, const intT &x) const { return x != zero; }
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
					res = max(res, value_type(abs(it.value())));
				}

			return std::move(res);
		}

		template <typename valT>
		bool is_zero(const SparseVector<valT> &v)
		{
			for (typename SparseVector<valT>::InnerIterator it(v); it; ++it)
			{
				if (it.value() != valT(0))
				{
					return false;
				}
			}

			return true;
		}

		template <typename valT>
		bool is_upper_triang(const SparseMatrix<valT> &A)
		{
			for (size_t k = 0; k < A.outerSize(); ++k)
				for (typename SparseMatrix<valT>::InnerIterator it(A, k); it; ++it)
				{
					if ((it.value() != 0) && (it.row() > it.col()))
						return false;
				}

			return true;
		}

		template <typename vecT>
		string dense_vec_to_str(const vecT &v)
		{
			stringstream res;

			if (!v.empty())
			{
				res << v[0];
				for (size_t i = 1; i < v.size(); ++i)
				{
					res << " " << v[i];
				}
			}

			return res.str();
		}

		// this functions transforms a sparce vector to a dense vector
		// and then prints it to a string
		template <typename valT>
		string sparse_vec_to_str(const SparseVector<valT> &v)
		{
			return dense_vec_to_str(conv<SparseVector<valT>, vector<valT>>::make(v));
		}
	} // namespace Utility
} // namespace ILPCuts
