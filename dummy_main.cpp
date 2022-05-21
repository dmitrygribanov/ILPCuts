// it's a sandbox
#include <iostream>

#include <eigen3/Eigen/SparseCore>
#include <NTL/ZZ_p.h>
#include <NTL/lzz_p.h>
#include <NTL/ZZ.h>

#include "linbox_module.h"
#include "ntl_module.h"
#include "flint_module.h"
#include "primes_downloader.h"

using namespace ILPCuts;
using namespace std;

#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/det.h>
#include <linbox/solutions/solve.h>

//#include <gmp++/gmp++.h>

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

    std::vector<int64_t> res = Utility::primes::download<int64_t>(10);
    for (auto p : res)
    {
        cout << p << " ";
    }
    cout << endl;

    // using Field1 = Givaro::Modular<Givaro::Integer>;
    // Field1 F1(101);

    // LinBox::SparseMatrix<Field1> A1(F1);
    // A1.read(std::cin);
    // Field1::Element d1;
    // LinBox::Method::Blackbox M;
    // det(d1, A1, M);

    // using Field2 = Givaro::QField<Givaro::Rational>;
    // Field2 F2;
    // LinBox::SparseMatrix<Field2> A2(F2);
    // Field2::Element d2;
    // det(d2, A2, M);

    Eigen::SparseMatrix<int64_t> A(2, 2);
    A.insert(0, 0) = 1;
    A.insert(0, 1) = 2;
    A.insert(1, 0) = 2;
    A.insert(1, 1) = 1;

    int64_t det = abs(Utility::NTL::det(A));

    Eigen::SparseMatrix<int64_t> H = Utility::NTL::compute_HNF(A, det);
    Eigen::SparseMatrix<int64_t> H1;
    Utility::flint::compute_HNF(A, H1);

    cout
        << H << endl;

    cout
        << H1 << endl;

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

    // ILPCuts::ConicCutsGenerator_adj<fmpzxx, fmpqxx, fmpzxx, fmpqxx> cuts_gen;

    // Eigen::SparseMatrix<fmpzxx> m(5, 5);

    // fmpzxx x, y;
    // fmpqxx q(1l, 2ul);
    // x = 7u;
    // y = x * x;
    // std::cout << x << "^2 = " << y << std::endl;

    return 0;
}