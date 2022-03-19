#pragma once

#include"conic_cuts_base.h"

namespace ILPCuts {

	template<typename intT, typename ratT, typename lintT, typename lratT>
	class ConicCutsGeneratorOpt 
		:public ConicCutsGenerator_Base<intT, ratT, lintT, lratT> {
	private:
		using base = typename ConicCutsGenerator_Base<intT, ratT, lintT, lratT>;
		using my_type = typename ConicCutsGeneratorOpt<intT, ratT, lintT, lratT>;

	public:
		using typename base::rational_type;
		using typename base::int_type;
		using typename base::lrational_type;
		using typename base::lint_type;

		using typename base::is_mat;
		using typename base::rs_mat;
		using typename base::lis_mat;
		using typename base::lrs_mat;

		using typename base::is_vec;
		using typename base::rs_vec;
		using typename base::lis_vec;
		using typename base::lrs_vec;

		using typename base::sparce_transform;
		using typename base::dir_coeffs;
		using typename base::process_func;
		using typename base::cone_traits;

		using base::tolerance;
		using base::threads;
		using base::out;
		using base::err;
		using base::err_tolerance;
		using base::cone_limit_per_direction;
		using base::direction_strategies;
		using base::_rat_to_int;
		using base::_lrat_to_lint;
		using base::rgen;
		using base::_to_single_int;
		using base::_to_single_rat;

		using base::base;

	private:		

		using perm_mat = SparseMatrix<int>;
		using LU_mat_value_type = lrational_type;
		using L_mat = TriangularView<SparseMatrix<LU_mat_value_type>, UpLoType::Lower>;
		using U_mat = TriangularView<SparseMatrix<LU_mat_value_type>, UpLoType::Upper>;

		struct LU_decomp_traits {
			L_mat L;
			U_mat U;
			perm_mat rowsP;
			perm_mat colsP;
		};

		struct cone_traits_spec :public base::cone_traits {

			using decomp_type = SparseLU<lrs_mat, COLAMDOrdering<int>>;
			std::unique_ptr<decomp_type> decomp;

		};

		std::unique_ptr<cone_traits> compute_cone_traits(const cone_traits* base_C,
			const sparce_transform& trans,
			std::vector<dir_coeffs> dir_coeffs_pull) override {

			auto new_cone = std::make_unique<cone_traits_spec>();
			new_cone->move_from_base_ptr(base::compute_cone_traits(base_C, trans, dir_coeffs_pull));

			try {
				new_cone->decomp = std::make_unique<cone_traits_spec::decomp_type>(new_cone->A.cast<lrational_type>());
			}
			catch (ToleranceErrorException exc) {
				out() << "Decomp. recompute tolerance error: " << exc.what() << '\n';
				err() << "Decomp. recompute tolerance error: " << exc.what() << '\n';
			}

			return std::move(new_cone);
		}

		std::unique_ptr<cone_traits_spec> create_cone_traits(const is_mat& A, const lint_type& det,
			const std::vector<rs_vec>& directions_pool) {
			auto res = std::make_unique<cone_traits_spec>();

			res->A = A;
			res->det = det;
			res->decomp = std::make_unique<cone_traits_spec::decomp_type>(A.cast<lrational_type>());

			for (size_t i = 0; i < directions_pool.size(); ++i) {
				auto t = Utility::system_solve(directions_pool[i].cast<lrational_type>().eval(), *(res->decomp) );
				res->dir_coeffs_pull.emplace_back(dir_coeffs{
					t.unaryExpr(_to_single_rat), i
					});
			}

			return std::move(res);
		}

		/*cone_traits construct_cone_traits(const is_mat & A) {
			using LU_solver = SparseLU<lint_type, COLAMDOrdering<int>>;

			cone_traits res;
			res.A = A;

			LU_solver solver(A);

			res.L = solver.matrixL();
			res.U = solver.matrixU();
			res.rowsP = solver.rowsPermutation();
			res.colsP = solver.colsPermutation();
			res.det = (*_lrat_to_lint)( abs(solver.determinant()) );

			return std::move(res);
		}*/

		void process_node(const cone_traits& base_C, process_func pf) override {

			const cone_traits_spec& C = static_cast<const cone_traits_spec&>(base_C);

			if (C.det % 2 == 0) {
				out() << "even step..." << endl;
				process_even_case(C, pf);
			}
			else {
				out() << "odd step..." << endl;
				process_odd_case(C, pf);
			}
		}

		is_vec even_system_solve(const cone_traits_spec& C) {
			size_t n = C.A.cols();

			is_vec res(n);

			bool odd_component_generated = false;
			while (!odd_component_generated) {

				//generation of a random rhs
				lrs_vec rhs = Utility::bernouly_vec_generate<int_type>(n, *rgen);

				/*cout << rhs << endl;
				cin.ignore();*/

				//solve the system
				lrs_vec x = std::move(Utility::system_solve(rhs, *C.decomp));

				//cout << "rhs : " << rhs << endl;
				//cout << "X : " << x << endl;

				//lrs_vec Z = C.A.cast<lrational_type>() * x - rhs;
				////cout <<"Z : " << Z << endl;
				//cout << "Z norm : " << Z.norm() << endl;

				/*cout << x << endl;
				cin.ignore();*/

				//create long integer solution
				lis_vec int_x = (x * C.det).unaryExpr(*_lrat_to_lint);

				lint_type max_2_pow = Utility::vec_max2pow(int_x);
				if ((static_cast<lint_type>(C.det) / max_2_pow) % 2 == 0) {
					int_x /= max_2_pow;
					//out() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					odd_component_generated = true;

					/*cout << int_x << endl;
					cin.ignore();*/

					//find odd solutions and convert them to 1 in result vector
					res = int_x.unaryExpr(std::ref(Utility::mod2_unaryop<lint_type>)).unaryExpr(_to_single_int);

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

		lis_vec odd_system_solve(const cone_traits_spec& C) {
			size_t n = C.A.cols();

			lis_vec res(n);

			bool nonzero_component_generated = false;
			while (!nonzero_component_generated) {

				//generation of a random rhs
				lrs_vec rhs = Utility::bernouly_vec_generate<int_type>(n, *rgen);

				//solve the system
				lrs_vec x = std::move(Utility::system_solve(rhs, *C.decomp));

				//create long integer solution
				res = (x * C.det).unaryExpr(*_lrat_to_lint);

				//transform modulo det
				res = res.unaryExpr(Utility::mod_unaryop<lint_type>(C.det));
				Utility::prune_zeroes_int(res);

				if (res.nonZeros() > 0) {
					nonzero_component_generated = true;
					cout << "res odd: " << res << endl;
				}
			}

			return std::move(res);
		}

		void process_even_case(const cone_traits_spec& C, process_func pf) {
			lis_vec x = even_system_solve(C).cast<lint_type>();
			for (typename lis_vec::InnerIterator it(x); it; ++it)
				if (it.value() > 0)
					pf(sparce_transform{
							static_cast<size_t>(it.index()), std::move(x), 2
						});
		}

		void process_odd_case(const cone_traits_spec& C, process_func pf) {
			lis_vec x = odd_system_solve(C);
			for (typename lis_vec::InnerIterator it(x); it; ++it)
				if (it.value() > 0) {
					lint_type lambda = Utility::mod_equation_solve<lint_type>(it.value(), C.det - gcd(C.det, it.value()), C.det);
					lis_vec t = (x * lambda).unaryExpr(Utility::mod_unaryop(C.det));

					pf(sparce_transform{
							static_cast<size_t>(it.index()), std::move(t), C.det
						});
				}
		}

	public:

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP, const lint_type& det) {

			auto directions_pool = std::move(base::generate_directions_pool(A));

			return std::move(base::generate_cuts(
				create_cone_traits(A, det, directions_pool),
				optLP
			));
		}

		std::vector<is_vec> generate_cuts(const is_mat& A, const rs_vec& optLP) override {

			auto start = std::chrono::steady_clock::now();

			lint_type det = (*_lrat_to_lint)( Utility::compute_abs_det(A.cast<lrational_type>().eval()) );

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			out() << "ConicCuts: determinant computation time: " << elapsed_seconds.count() << "s\n";

			return std::move(generate_cuts(A, optLP, det));
		}
	};
}
