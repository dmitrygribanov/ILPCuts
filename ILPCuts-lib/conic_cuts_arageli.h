#pragma once

#include<arageli.hpp>
#include<vector>

template<typename intT, typename rationalT>
class ConicCutsGenerator_Arageli {
public:

	typedef rationalT rational_type;
	typedef intT int_type;
	typedef Arageli::vector<rational_type> rational_vec;
	typedef Arageli::matrix<rational_type> rational_mat;
	typedef Arageli::vector<int_type> int_vec;
	typedef Arageli::matrix<int_type> int_mat;

	/*class ChoiceStrategy {
	protected:
		std::vector<size_t> _base_idxs;
		const int_mat* _facets_mat;
		const int_vec* _cT;

	public:
		ChoiceStrategy() 
			:_facets_mat(nullptr), _cT(nullptr) 
		{}

		virtual void init(const int_vec& cT, const int_mat& facets_mat) {
			_facets_mat = &facets_mat;
			_cT = &cT;

			size_t n = _facets_mat->nrows();
			_base_idxs.resize(n);
			for (size_t i = 0; i < n; ++i)
				_base_idxs[i] = i;

#ifdef _DEBUG
			{
				int_mat A = _facets_mat->copy_cols(_base_idxs);
				assert(Arageli::det(A) != 0, "First n rows of facets_mat must be LI.");
			}
#endif 
		}

		virtual bool next_base() { return false; }
		 
		const std::vector<size_t>& get_base_idxs()const {
			return _base_idxs;
		}
	};*/

	rational_mat make_cuts(const rational_vec& vertex_LP,
		const int_vec& cT,
		int_mat incident_facets,
		std::unique_ptr<ChoiceStrategy> strategy = std::make_unique<ChoiceStrategy>() )
	{
		size_t n = vertex_LP.size();
#ifdef _DEBUG
 		assert(incident_facets.nrows() >= n, 
			"Facets number must be >= dim(vertex_LP)");
		assert(incident_facets.ncols() == n + 1, 
			"Rows of incident_facets must have (dim(vertex_LP)+1) components");
#endif

		int_mat A(incident_facets.nrows(), n, 0);
		int_vec b(incident_facets.nrows(), 0);

		for (size_t i = 0; i < incident_facets.nrows(); ++i) {
			for (size_t j = 0; j < n; ++j)
				A(i, j) = incident_facets(i, j);

			b[i] = incident_facets(i, n);
		}

		A.transpose();

#ifdef _DEBUG
		assert(Arageli::rank(A) == n, 
			"rank(incident_facets) must be = dim(vertex_LP)");
		assert(vertex_LP * A == b,
			"(A * vertex_LP) must be equal to b");
#endif

		/*Arageli::output_aligned(std::cout, A, "| | ", " | |", " ");
		Arageli::output_aligned(std::cout, int_mat(b), "| | ", " | |", " ");*/
		_normalize(A, b);
		/*Arageli::output_aligned(std::cout, A, "| | ", " | |", " ");
		Arageli::output_aligned(std::cout, int_mat(b), "| | ", " | |", " ");*/

		int_mat all_cuts(0,n+1,0);
		strategy->init(cT,A);
		do {
			int_mat A_base = A.copy_cols(strategy->get_base_idxs());
			//int_vec b_base = b.copy_subvector(strategy->get_base_idxs());

			//Arageli::output_aligned(std::cout, A_base, "| | ", " | |", " ");
			//Arageli::output_aligned(std::cout, b_base, "| | ", " | |", " ");

			int_mat base_cuts = _make_cuts_from_square_mat(vertex_LP, cT, A_base );
			
			all_cuts.insert_rows(all_cuts.nrows(), base_cuts.nrows(), base_cuts.begin() );

			/*Arageli::output_aligned(std::cout, base_cuts, "| | ", " | |", " ");
			Arageli::output_aligned(std::cout, all_cuts, "| | ", " | |", " ");*/

		} while ( strategy->next_base() );

		return all_cuts;
	}

private:

	void _normalize(int_mat & A, int_vec& b) {
		rational_mat res, transform;
		std::vector<size_t> basis;
		Arageli::rref_gauss_field(A, res, transform, basis);

		for (size_t i = 0; i < basis.size(); ++i) {
			A.swap_cols(i, basis[i]);
			b.swap_els(i, basis[i]);
		}
	}

	int_mat _build_unimodular_cone(const int_vec& cT, const int_mat& A) 
	{
		size_t n = A.ncols();
		int_mat S, P, Q;
		int_type det;
		size_t rank;
		Arageli::smith(A, S, P, Q, rank, det);
		const int_type& Snn = S(n - 1, n - 1);


		//Arageli::output_aligned(std::cout, S, "| | ", " | |", " ");
		//Arageli::output_aligned(std::cout, P, "| | ", " | |", " ");
		//Arageli::output_aligned(std::cout, Q, "| | ", " | |", " ");

		if( Snn == 1 ) {
			return A;
		}
		else {
			int_vec Qn = Q.copy_col(n - 1);

			//std::cout << Qn << '\n';

			std::vector<int_mat> cone_candidates;

			if (Snn % 2 == 0) {
				//int_vec t_int = ((Snn / 2) * Qn) % Snn;
				rational_vec t = rational_vec(Qn % 2) / 2;
				/*for (size_t i = 0; i < n; ++i)
					t[i] = (Qn[i] % 2) / 2;*/
				//std::cout << t << '\n';

				for (size_t i = 0; i < n; ++i) {
					if (t[i] != 0) {
						int_mat cone_candt = A;
						cone_candt.assign_col(i, A * t);
						cone_candidates.push_back(cone_candt);
					}
				}
			} else {
				for (size_t i = 0; i < n; ++i) {
					if (Qn[i] != 0) {
						int_mat cone_candt = A;

						int_type lambda = _mod_equation_solve(Qn[i], Snn - Arageli::gcd(Qn[i], Snn), Snn);
						int_vec t_int = (lambda * Qn) % Snn;
						rational_vec t = rational_vec(t_int) / Snn;

						cone_candt.assign_col(i, A * t);
						cone_candidates.push_back(cone_candt);
					}
				}
			}

			for (auto& cone_candt : cone_candidates) {
				rational_vec x = Arageli::solve_linsys_field(rational_mat(cone_candt), rational_vec(cT));
				if (x >= 0) {
					return _build_unimodular_cone(cT, cone_candt);
				}
			}
		}
	}

	int_type _mod_equation_solve(int_type a, int_type b, int_type mod) const {
#ifdef _DEBUG
		assert(mod >= 0 && a >= 1 && b >= 0 && b < mod);
#endif
		int_type x = 0, y = 0;
		int_type gcd_a_mod = Arageli::euclid_bezout(a, mod, x, y);

		if (x < 0) x += mod;

#ifdef _DEBUG
		assert(b % gcd_a_mod == 0, "No solutions in _mod_equation_solve.");
#endif
		x *= (b / gcd_a_mod);

		return x;
	}

	int_mat _make_cuts_from_square_mat(const rational_vec& vertex_LP,
		const int_vec& cT,
		int_mat& A) 
	{
#ifdef _DEBUG
		{
			assert(A.ncols() == A.nrows(), "A must be square.");
			assert(A.nrows() == cT.size() && cT.size() == vertex_LP.size(),
				"Dimensions of A, cT, vertex_LP must be the same.");

			rational_vec x;
			Arageli::solve_linsys_result res = 
				Arageli::solve_linsys_field(rational_mat(A), rational_vec(cT), x);

			assert(res == Arageli::solve_linsys_result::SLR_UNIQUE, "A must be non-singular.");
			assert(x >= 0, "cT must be in cone(A.transpose()).");
		}
#endif

		int_mat unimodular_mat = _build_unimodular_cone(cT, A );
		unimodular_mat.transpose();

		int_vec b = Arageli::floor(unimodular_mat * vertex_LP);

		unimodular_mat.insert_col(unimodular_mat.ncols(), b);

		return unimodular_mat;
	}
};
