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

namespace ILPCuts {

	using namespace Eigen;

	class ToleranceErrorException : public std::exception {
	public:
		ToleranceErrorException()
			:std::exception() {}

		ToleranceErrorException(std::string msg)
			:std::exception(msg.c_str()) {}
	};

	namespace Utility {
		template<typename intT, typename ratT>
		intT rat_to_int(const ratT& rat_x, const ratT& err_tolerance) {
			//std::cout << rat_x << std::endl;
			intT int_x = static_cast<intT>(round(rat_x));
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

		template<typename intT>
		intT mod2_unaryop(const intT& x) {
			return abs(x % 2);
		}

		template<typename intT>
		class mod_unaryop {
		private:
			const intT& _mod;

		public:
			mod_unaryop(const intT& mod)
				:_mod(mod) {}

			intT operator()(const intT& x) const {
				intT res = x % _mod;
				if (res < 0) res += _mod;
				return res;
			}
		};

		template<typename ratT>
		class Round {
		public:
			ratT operator()(const ratT& rat_x) const {
				return std::move(round(rat_x));
			}
		};

		template<typename ratT>
		class Floor {
		public:
			ratT operator()(const ratT& rat_x) const {
				return std::move(floor(rat_x));
			}
		};

		template<typename intT, typename ratT>
		class RatToInt_Unary {
		private:
			ratT _err_tolerance;

		public:
			RatToInt_Unary(const ratT& err_tolerance)
				:_err_tolerance(err_tolerance) {}

			intT operator()(const ratT& rat_x) const {
				return rat_to_int<intT, ratT>(rat_x, _err_tolerance);
			}
		};

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

		template<typename ratMatrixT, typename solverT>
		ratMatrixT compute_inverse(const ratMatrixT& A, solverT& solver);

		template<typename rationalT, typename solverT>
		SparseMatrix<rationalT> compute_inverse(const SparseMatrix<rationalT>& A, solverT& solver);

		/*template<typename rationalT, typename ratMatrixT>
		static std::tuple<ratMatrixT, rationalT> compute_adjugate_and_absdet(const ratMatrixT& A);*/

		template<typename rationalT>
		std::tuple<SparseMatrix<rationalT>, rationalT> compute_abs_adjugate_and_det(const SparseMatrix<rationalT>& A);

		/*template<typename ratVectorT, typename ratMatrixT>
		static ratMatrixT adjust_normal_cone(const ratMatrixT& A, const ratVectorT& c);

		template<typename rationalT>
		static SparseMatrix<rationalT> adjust_normal_cone(const SparseMatrix<rationalT>& A, const SparseVector<rationalT>& c, const rationalT& tolerance);*/

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

		template<typename intT>
		intT mod_equation_solve(const intT& a, const intT& b, const intT& mod);

		template<>
		mpz_class mod_equation_solve(const mpz_class& a, const mpz_class& b, const mpz_class& mod) {
#ifdef _DEBUG
			assert(mod > 0 && a >= 1 && b >= 0 && b < mod);
#endif

			mpz_class x = 0, y = 0;
			mpz_class gcd_a_mod = 0; // = Internal::euclid_bezout(a, mod, x, y);

			mpz_gcdext(gcd_a_mod.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

			if (x < 0) x += mod;

#ifdef _DEBUG
			assert(x > 0);
			assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
			x *= (b / gcd_a_mod);

			return x;
		}
	};

	template<typename rationalT, typename solverT>
	SparseMatrix<rationalT> Utility::compute_inverse(const SparseMatrix<rationalT>& A, solverT& solver) {
		using rat_mat = SparseMatrix<rationalT>;

		size_t n = A.cols();

		rat_mat I(n, n);
		I.reserve(VectorXi::Constant(n, 2));
		for (size_t i = 0; i < n; ++i)
			I.insert(i, i) = 1;

		rat_mat res(n, n);
		for (size_t i = 0; i < n; ++i) {
			res.col(i) = solver.solve(I.col(i));
		}

		return std::move(res);
	}

	template<typename rationalT>
	std::tuple<SparseMatrix<rationalT>, rationalT> Utility::compute_abs_adjugate_and_det(const SparseMatrix<rationalT>& A) {

		using rat_type = rationalT;
		using rat_mat = SparseMatrix<rat_type>;

		//SparseLU<rs_mat> decomp_LU(A.cast<rational_type>());
		SparseLU<rat_mat, COLAMDOrdering<int>> decomp_LU(A);

		rat_type det = abs(decomp_LU.determinant());
		rat_mat adjA = compute_inverse<rat_type>(A, decomp_LU) * det;

		/*size_t n = A.cols();
		rat_mat I(n, n);
		I.reserve(VectorXi::Constant(n, 2));
		for (size_t i = 0; i < n; ++i)
			I.insert(i, i) = det;*/

		return std::make_tuple( std::move(adjA), std::move(det) );
	}

	template<typename intT, typename rationalT>
	class ConicCutsGenerator_Base {
	public:
		using rational_type = rationalT;
		using int_type = intT;

		using is_mat = SparseMatrix<int_type>;
		using rs_mat = SparseMatrix<rational_type>;

		using is_vec = SparseVector<int_type>;
		using rs_vec = SparseVector<rational_type>;

	private:

		std::ostream* _out;
		std::ostream* _err;

	protected:
		using RatToInt_Unary = typename Utility::RatToInt_Unary<int_type, rational_type>;
		std::unique_ptr<RatToInt_Unary> _rat_to_int;

	public:

		class DirectionsStrategyBase {
		public:
			virtual std::vector<rs_vec> generate_directons(const is_mat& A) = 0;
		};

		rational_type tolerance;
		rational_type err_tolerance;
		size_t cone_limit_per_direction;
		size_t threads;
		std::vector<std::unique_ptr<DirectionsStrategyBase> > direction_strategies;

		ConicCutsGenerator_Base(std::ostream* outs = &std::cout, std::ostream* errs = &std::cout)
			:_out(outs), _err(errs), cone_limit_per_direction(0),
			tolerance(0), err_tolerance(0), threads(0), _rat_to_int() 
		{}

		void Init(const rational_type& a_tolerance = 0.000001,
			const rational_type& a_err_tolerance = 0.001,
			size_t a_cone_limit_per_direction = 1,
			size_t a_threads = 8) {
			tolerance = a_tolerance;
			err_tolerance = a_err_tolerance;
			cone_limit_per_direction = a_cone_limit_per_direction;
			direction_strategies.push_back(std::make_unique<SomeDirection>());
			threads = a_threads;

			_rat_to_int = std::make_unique<RatToInt_Unary>(RatToInt_Unary(err_tolerance));
		}

		std::ostream& out() { return *_out; }

		std::ostream& err() { return *_err; }

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

		is_mat adjust_normal_cone(const is_mat& A, const is_vec& c) {
			is_mat resA = A;

			rs_mat adjA;
			rational_type det;

			std::tie(adjA, det) = Utility::compute_abs_adjugate_and_det<rational_type>(A.cast<rational_type>());

			rs_vec t = adjA * c.cast<rational_type>();
			/*cout << c << endl << endl;
			cout << t << endl;
			cin.ignore();*/

			for (typename rs_vec::InnerIterator it(t); it; ++it) {
				// it.value() < tolerance previouse code. working fast but its must be not correct
				if ( it.value() < 0 && abs(it.value()) > tolerance ) {
					/*cout << it.value() << endl;
					cin.ignore();*/
					resA.col(it.index()) = -resA.col(it.index());
				}
			}

			return std::move(resA);
		}

	};
}
