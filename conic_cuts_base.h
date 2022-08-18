#pragma once

#include "base_module.h"

// all optional modules are added at the moment
// to simplify experiments
#include "Eigen_solving_utilities.h"

namespace ILPCuts
{

	template <typename intT, typename rationalT>
	class cuts_generator_base
	{
	public:
		using rational_type = rationalT;
		using int_type = intT;

		using is_mat = SparseMatrix<int_type>;
		using rs_mat = SparseMatrix<rational_type>;

		using is_vec = SparseVector<int_type>;
		using rs_vec = SparseVector<rational_type>;

		using Strategies = CutsGenStrategies<int_type, rational_type>;

		struct cut_info
		{
			is_vec ray;
			int_type rhs;
			rational_type frac_rhs;

			bool operator<(const cut_info &arg)
			{
				bool is_rays_less = Utility::SparseVecLess<int_type>()(ray, arg.ray);
				if (is_rays_less)
					return true;

				bool is_rays_equal = Utility::SparseVecEquiv<int_type>()(ray, arg.ray);
				if (is_rays_equal)
				{
					return rhs < arg.rhs;
				}

				return false;
			}

			bool operator==(const cut_info &arg)
			{
				return Utility::SparseVecEquiv<int_type>()(ray, arg.ray) && rhs == arg.rhs;
			}
		};

	protected:
		using RatToInt_Unary = Utility::RatToInt_Unary<int_type, rational_type>;

		RatToInt_Unary _rat_to_int;

		is_mat _AT_copy;
		rs_vec _optLP_copy;

		std::vector<cut_info> _resulting_cuts;

		void copy_original_problem_info(const is_mat &A, const rs_vec &optLP)
		{
			_AT_copy = A;
			_optLP_copy = optLP;
		}

		vector<cut_info> generate_cuts_from_row_matrix(const is_mat &cuts_mat, const rs_vec &optLP)
		{
			vector<cut_info> res;

			rs_vec rat_b = cuts_mat.template cast<rational_type>() * optLP;
			is_vec rounded_b = rat_b.unaryExpr(Utility::Round_Unary<rational_type>()).unaryExpr(_rat_to_int);
			is_vec floor_b = rat_b.unaryExpr(Utility::Floor_Unary<rational_type>()).unaryExpr(_rat_to_int);

			rs_vec diff_b = (rat_b - rounded_b.template cast<rational_type>()).cwiseAbs();

			typename rs_vec::InnerIterator diff_it(diff_b);
			while (diff_it)
			{
				if (diff_it.value() > rational_type(tolerance))
				{
					is_vec ray = cuts_mat.row(diff_it.index());
					int_type rhs = floor_b.coeff(diff_it.index());
					rational_type frac_rhs = rat_b.coeff(diff_it.index());

					res.emplace_back(cut_info{ray, rhs, frac_rhs});
				}
				++diff_it;
			}

			return res;
		}

		std::vector<rs_vec> generate_directions_pool(const is_mat &A)
		{
			std::vector<rs_vec> directions_pool;

			for (const auto &strategy : direction_strategies)
			{
				std::vector<rs_vec> strategy_directions = std::move(strategy->generate_directons(A));
				directions_pool.insert(directions_pool.end(),
									   std::make_move_iterator(strategy_directions.begin()),
									   std::make_move_iterator(strategy_directions.end()));
			}

			return std::move(directions_pool);
		}

	public:
		logger &log;

		double tolerance = 0.000001;
		double err_tolerance = 0.01;
		size_t threads = 8;
		std::mt19937 *rgen = nullptr;
		size_t cone_limit_per_direction = 1;

		std::vector<std::unique_ptr<typename Strategies::DirectionsStrategyBase>> direction_strategies;

		cuts_generator_base(logger &a_log)
			: log(a_log), _rat_to_int(err_tolerance)
		{
		}

		logger::stream &out()
		{
			log.out() << "ILPCuts-Generator: ";
			return log.out();
		}

		const std::vector<cut_info> &get_cuts() const { return _resulting_cuts; }

		is_mat adjust_normal_cone(is_mat A, const rs_vec &c)
		{
			rs_vec t = Utility::system_solve(A.template cast<rational_type>().eval(), c);

			for (typename rs_vec::InnerIterator it(t); it; ++it)
			{
				// it.value() < tolerance previouse code. working fast but its must be not correct
				if (it.value() < 0 && abs(it.value()) > tolerance)
				{
					A.col(it.index()) = -A.col(it.index());
				}
			}

			return std::move(A);
		}

		// this function takes a collection of resulting cuts
		// and prints them to a string
		// each cut is separtated by \n
		// the last cut is also finished by \n
		void export_cuts(ofstream &out_file)
		{
			for (size_t i = 0; i < _resulting_cuts.size(); ++i)
			{
				out_file << Utility::vec_to_str_dense(_resulting_cuts[i].ray) << " "
						 << _resulting_cuts[i].rhs << "\n";
			}
		}

		void clear_cuts()
		{
			_resulting_cuts.clear();
			_AT_copy.resize(0, 0);
			_optLP_copy.resize(0);
		}

		void check_cuts_validity() const
		{
			size_t n = _AT_copy.cols();
			assert(n > 0);

			for (size_t i = 0; i < _resulting_cuts.size(); ++i)
			{
				const auto &cut = _resulting_cuts[i];

				rs_vec t = Utility::system_solve(_AT_copy.template cast<rational_type>().eval(), cut.ray.template cast<rational_type>().eval());

				for (typename rs_vec::InnerIterator it(t); it; ++it)
				{
					if (it.value() < 0 && abs(it.value()) > tolerance)
					{
						stringstream msg;
						msg << "ILPCuts error: cut(" << i << ") is invalid: negative solution(" << it.value() << ")...\n";
						throw Exception(msg.str());
					}
				}

				if (abs(cut.ray.template cast<rational_type>().dot(_optLP_copy) - cut.frac_rhs) > err_tolerance)
				{
					stringstream msg;
					msg << "ILPCuts error: cut(" << i << ") is invalid: invalid fractional rhs...\n";
					throw Exception(msg.str());
				}

				if (floor(cut.frac_rhs) != cut.rhs)
				{
					stringstream msg;
					msg << "ILPCuts error: cut(" << i << ") is invalid: invalid floor rounding...\n";
					throw Exception(msg.str());
				}
			}
		}

		virtual void generate_cuts(const is_mat &A, const rs_vec &optLP) = 0;
	};
} // namespace ILPCuts