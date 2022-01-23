#pragma once

#include<numeric>

namespace ILPCuts {

	namespace Internal {
		/*template <typename T>
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
		}*/
	}

	/*template<typename intT>
	intT gcd(const intT& a, const intT& b) {
		return std::gcd(a, b);
	}*/

	/*template <typename T>
	T euclid(T a, T b)
	{
		while ((a %= b) != 0)
		{
			std::swap(a, b);
		}
		return b;
	}

	template <typename T>
	inline T gcd
	(
		const T& a,
		const T& b
	)
	{
		return euclid(a, b);
	}*/

}
