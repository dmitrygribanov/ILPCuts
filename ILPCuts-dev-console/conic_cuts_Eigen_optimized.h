#pragma once

#include"base_module.h"

namespace ILPCuts {

	using namespace Eigen;

	template<typename intT, typename ratT, typename lintT, typename lratT>
	class ConicCutsGeneratorOpt :public ConicCutsGenerator_Base<intT, ratT> {
	private:
		using base = typename ConicCutsGenerator_Base<intT, rationalT>;

	public:
		using typename base::rational_type;
		using typename base::int_type;

		using typename base::is_mat;
		using typename base::rs_mat;

		using typename base::is_vec;
		using typename base::rs_vec;

		using lint_type = lintT;
		using lrational_type = lratT;

		using lis_mat = SparseMatrix<lint_type>;
		using lrs_mat = SparseMatrix<lrational_type>;

		using lis_vec = SparseVector<lint_type>;
		using lrs_vec = SparseVector<lrational_type>;

		using base::tolerance;
		using base::out;
		using base::err;
		using base::err_tolerance;
		using base::cone_limit_per_direction;
		using base::direction_strategies;

		using base::base;
	};
}
