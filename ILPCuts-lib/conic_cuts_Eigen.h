#pragma once

#include<Eigen\SparseCore>
#include<Eigen\Dense>

#include<vector>
#include<numeric>

namespace ILPCuts {

	using namespace Eigen;

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
		int_type vec_gcd(const int_vec& vec) const {
			int_type res = std::gcd(vec[0], vec[1]);

			for (size_t i = 2; i < vec.size(); ++i) {
				res = std::gcd(res, vec[i]);
			}

			return res;
		}

		int_type vec_max2pow(const int_vec& vec) const {
			int_type res = 1;

			while ( (vec.array().unaryExpr([](const int_type & x) { return x % 2; }) 
				== 0).all()  ) {
				res *= 2;
			}

			return res;
		}

		int_type is_even(const int_type& x) {

		}

	public:

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
				return rat_mat(_c);
			}
		};

		void partial_uni_decomp
		(std::vector<int_mat> & result, const int_mat& A, const rat_mat& directions_pool,
			const int_mat& adjA, const int_type& det,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {

			if (det == 1)
				result.push_back(A);
			else {
				size_t n = A.cols();

				if (det % 2 == 0) {
					int_vec u;
					for (size_t i = 0; i < n; ++i) {
						u = adjA.col(i);
						int_type adj_col_2pow = vec_max2pow(u);
						if (((int_type)det / adj_col_2pow) % 2 == 0) {
							u /= adj_col_2pow;
							break;
						}
					}

					rat_vec t(u.size(), 0);
					for (size_t i = 0; i < t.size(); ++i) {
						if (u[i] % 2 != 0)
							t[i] = 0.5;
					}

					for (size_t i = 0; i < n; ++i) {
						int_mat newA = A;
						newA.col(i) = (A.cast<rational_type>() * t).array().round().cast<int_type>();
					}
					
				} else {

				}
			}
		}

		std::vector<int_vec> generate_cuts(const int_mat& A, const int_mat& adjA, const rat_vec& optLP,
			DirectionsStrategyBase & directions_strategy,
			const int_type& det,
			const rational_type& tolerance = 0.00001,
			const rational_type& err_tolerance = 0.0001) {
			
			std::vector<int_vec> res_cuts;

			rat_mat directions_pool = directions_strategy.generate_directons();

			std::vector<int_mat> uni_mat_pool;
			
			partial_uni_decomp(uni_mat_pool, A, directions_pool,
				adjA, det, tolerance, err_tolerance);
			
			for (const int_mat& uni_mat : uni_mat_pool) {
				rat_vec rat_b = uni_mat.cast<rational_type>() * optLP;
				int_vec int_b = rat_b.array().floor().cast<int_type>();

				rat_vec diff_b = (rat_b - int_b.cast<rational_type>()).array().abs();

				for (size_t i = 0; i < diff_b.size(); ++i) {
					if (diff_b[i] > tolerance) {
						int_vec cut = uni_mat.row(i);
						cut << int_b[i];
						res_cuts.push_back(cut);
					}
				}
			}

			return std::move(res_cuts);
		}
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
