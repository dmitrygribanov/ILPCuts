#pragma once

#include <vector>
#include <random>
#include <atomic>
#include <string>
#include <sstream>
#include <future>
#include <thread>
#include <functional>

#include "numeric_module.h"
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
	namespace Utility
	{
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
