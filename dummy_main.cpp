// it's a sandbox
#include <iostream>

#include <eigen3/Eigen/SparseCore>
#include <NTL/ZZ_p.h>
#include <NTL/lzz_p.h>
#include <NTL/ZZ.h>

#include "base_module.h"

using namespace ILPCuts;
using namespace std;

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/det.h>

// #include "conic_cuts_Eigen.h"

int main(int argc, char **argv)
{
    // Eigen::SparseMatrix<NTL::zz_p> m(2, 2);
    // using Tripl = Eigen::Triplet<NTL::zz_p>;
    // Tripl item(0, 0, NTL::zz_p(5l));
    // std::vector<Tripl> triplets;
    // triplets.emplace_back(Tripl(0, 0, NTL::zz_p(5l)));
    // triplets.back().value().init(7l);
    // triplets.push_back(Tripl(0, 0, NTL::zz_p(2l)));
    // triplets.back().value().init(7l);

    // m.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<Utility::primes::ptype> res = Utility::primes::download(10);
    for (auto p : res)
    {
        cout << p << " ";
    }
    cout << endl;

    using Field = Givaro::Modular<int64_t>;
    Field F(101);

    LinBox::SparseMatrix<Field> A(F);
    A.read(std::cin);
    Field::Element d;
    LinBox::Method::Blackbox M;
    det(d, A, M);

    // nmodxx_ctx n(7ul);
    // nmodxx x = nmodxx::make_nored(3, n);
    // using Tripl = Eigen::Triplet<nmodxx>;
    // Tripl item(0, 0, nmodxx::make_nored(0, n));
    // std::vector<Tripl> triplets;
    // triplets.emplace_back(Tripl(0, 0, nmodxx::make_nored(0, n)));
    // triplets.push_back(Tripl(1, 1, nmodxx::make_nored(1, n)));
    // Eigen::SparseMatrix<nmodxx> m(2, 2);
    // m.setFromTriplets(triplets.begin(), triplets.end());
    // m.coeff(0, 0) = nmodxx::make_nored(0, n);
    // m.coeff(0, 1) = nmodxx::make_nored(1, n);
    // m.coeff(1, 0) = nmodxx::make_nored(2, n);
    // m.coeff(1, 1) = nmodxx::make_nored(3, n);

    // fmpqxx x, y, z;
    // z = x / y;

    // fmpzxx i;

    // fmpqxx k = (x * i).evaluate();
    // double dd;
    // x < fmpqxx(dd);

    // ILPCuts::ConicCutsGenerator<fmpzxx, fmpqxx, fmpzxx, fmpqxx> cuts_gen;

    // Eigen::SparseMatrix<fmpzxx> m(5, 5);

    // fmpzxx x, y;
    // fmpqxx q(1l, 2ul);
    // x = 7u;
    // y = x * x;
    // std::cout << x << "^2 = " << y << std::endl;

    return 0;
}