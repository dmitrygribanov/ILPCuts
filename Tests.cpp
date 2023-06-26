#include "wiedemann_alg.h"

#include <eigen3/Eigen/LU>

using namespace ILPCuts;
using namespace std;

namespace ILPCuts
{
    using namespace NTL;

    logger log("tests_log.txt");

    void seq_min_poly_unit_test()
    {
        log.out() << "Sequence minimal polynomial test..." << endl;

        zz_p::init(5);
        log.out() << "modulus: " << zz_p::modulus() << endl;

        stringstream str;
        zz_pX fX;
        str << "[3 0 4 2 3 0]";
        str >> fX;

        log.out() << "sequence: " << str.str() << endl;

        zz_pX mX = ModuleNTL::seq_min_poly(fX);
        log.out() << "minimal polynomial: " << mX << endl;

        str << "[2 2 1]";
        zz_pX ansX;
        str >> ansX;
        log.out() << "correct answer is: " << ansX << endl;

        if (mX == ansX)
            log.out() << "---DONE---\n"
                      << endl;
        else
            log.out() << "---FAILED---\n"
                      << endl;
    }

    void pade_approx_unit_test()
    {
        log.out() << "Pade's approximation test..." << endl;

        zz_p::init(5);
        log.out() << "modulus: " << zz_p::modulus() << endl;

        stringstream str;
        zz_pX fX;
        str << "[1 2 3 4]";
        str >> fX;

        log.out() << "sequence(x) = " << str.str() << "\n"
                  << "(k, n-k) = (2,2)" << endl;

        zz_pX sX, tX;
        tie(sX, tX) = ModuleNTL::pade_approx(fX, 4, 2);

        log.out() << "result is (s/t):\n"
                  << "s(x) = " << sX << "\n"
                  << "t(x) = " << tX << endl;

        zz_pX ans_sX(INIT_MONO, 0, 1);
        zz_pX ans_tX;

        str << "[1 3 1]";
        str >> ans_tX;

        log.out() << "correct answer is (s/t):\n"
                  << "s(x) = " << ans_sX << "\n"
                  << "t(x) = " << ans_tX << endl;

        if (sX == ans_sX && tX == ans_tX)
            log.out() << "---DONE---\n"
                      << endl;
        else
            log.out() << "---FAILED---\n"
                      << endl;
    }

    void vector_min_poly_unit_test()
    {
        stringstream str;
        log.out() << "Vector minimal polynomial test..." << endl;

        zz_p::init(5);
        log.out() << "modulus: " << zz_p::modulus() << endl;

        MatrixXi A_dense(3, 3);
        A_dense << 1, 4, 4,
            4, 0, 3,
            1, 2, 4;

        log.out() << "matrix A:\n"
                  << A_dense << endl;

        SparseMatrix<int64_t> A_sparse(A_dense.cast<int64_t>().sparseView());
        ModuleNTL::Oracle_zzp<int64_t> O(A_sparse);

        ModuleNTL::WiedemannsAlg wd_alg(log);

        str << "[3 1 2]";
        vec_zz_p b;
        str >> b;

        // std::copy(b.begin(), b.end(), std::ostream_iterator<char>(log.out(), " "));
        log.out() << "vector b: "
                  << b << endl;

        zz_pX mX = wd_alg.vec_min_poly(O, b);
        log.out() << "minimal polynomial: " << mX << endl;

        str << "[1 3 0 1]";
        zz_pX ans_mX;
        str >> ans_mX;

        if (mX == ans_mX)
            log.out() << "---DONE---\n"
                      << endl;
        else
            log.out() << "---FAILED---\n"
                      << endl;
    }

    void Zp_kernel_unit_test(const MatrixXi &A_dense, int modulus)
    {
        log.out() << "Zp kernel test..." << endl;

        log.out() << "matrix A:\n"
                  << A_dense << endl;

        FullPivLU<MatrixXd> LU(A_dense.cast<double>());
        int det = abs(round(LU.determinant()));

        zz_p::init(modulus);
        log.out() << "det(A) = " << det << endl;
        log.out() << "modulus = " << modulus << endl;

        SparseMatrix<int64_t> A_sparse(A_dense.cast<int64_t>().sparseView());
        ModuleNTL::Oracle_zzp<int64_t> O(A_sparse);

        ModuleNTL::WiedemannsAlg wd_alg(log);

        vec_zz_p dense_x = wd_alg.solve_kernel_v1(O);
        log.out() << "solution \"x\" is: " << dense_x << endl;

        auto sparse_x = conv<vec_zz_p, SparseVector<int64_t>>::make(dense_x);
        if (Utility::is_zero((A_sparse * sparse_x).unaryExpr(Utility::mod_unaryop(modulus)).eval()))
        {
            log.out() << "A x == 0\n"
                      << "---DONE---\n"
                      << endl;
        }
        else
        {
            log.out() << "A x != 0\n"
                      << "---FAILED---\n"
                      << endl;
        }
    }

} // namespace ILPCuts

// void unit_test_I() {
//
//     cout << "---------------unit_test_I---------------------" << endl;
//     using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<long long int, double>;
//     CutsGen gen;
//
//     CutsGen::int_mat        A(4,4);
//     CutsGen::rat_vec      optLP(4);
//     CutsGen::rat_vec          c(4);
//     CutsGen::SpecifiedDirection gd(c);
//
//     A << 1,   0,   0, 0,
//          1,   1,   0, 0,
//         10, -10,   1, 0,
//         24,  30, -30, 1;
//
//     CutsGen::int_vec t(4);
//     t << 1, 1, 1, 1;
//
//     c = (A * t).cast<double>();
//
//     optLP << 1, 2, 4, -5.5;
//
//     std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);
//
//
//     assert(cuts.empty());
//
//     cout << "-----------------------------------------------" << endl;
// }
//
// void unit_test_S2() {
//     cout << "---------------unit_test_S2---------------------" << endl;
//     using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<long long int, double>;
//     CutsGen gen;
//
//     CutsGen::int_mat A(3, 3);
//     std::vector<CutsGen::int_vec> true_cuts;
//     CutsGen::rat_vec optLP(3);
//     CutsGen::rat_vec c(3);
//     CutsGen::SpecifiedDirection gd(c);
//
//     A << 1, 0, 0,
//          1, 1, 0,
//          0, 0, 2;
//
//     c << 1, 1, 1;
//
//     optLP << 1, 2, 5.5;
//
//     std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);
//
//     true_cuts.push_back(CutsGen::int_vec(4));
//
//     true_cuts[0] << 0, 0, 1, 5;
//
//
//     assert( cuts == true_cuts );
//
//     cout << "-----------------------------------------------" << endl;
// }
//
// void unit_test_S3() {
//     cout << "---------------unit_test_S3---------------------" << endl;
//     using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<long long int, double>;
//     CutsGen gen;
//
//     CutsGen::int_mat A(3, 3);
//     std::vector<CutsGen::int_vec> true_cuts;
//     CutsGen::rat_vec optLP(3);
//     CutsGen::rat_vec c(3);
//     CutsGen::SpecifiedDirection gd(c);
//
//     A << 1, 0, 0,
//         1, 1, 0,
//         0, 0, 3;
//
//     c << 1, 1, 1;
//
//     optLP << 1, 2, (float)(11) / 3;
//
//     std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);
//
//     true_cuts.push_back(CutsGen::int_vec(4));
//
//     true_cuts[0] << 0, 0, 1, 3;
//
//     assert(cuts == true_cuts);
//
//     cout << "-----------------------------------------------" << endl;
// }

// void unit_test_POW2() {
//     cout << "---------------unit_test_POW2---------------------" << endl;
//     using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<long long int, double>;
//     CutsGen gen;
//
//     CutsGen::int_mat A(3, 3);
//     std::vector<CutsGen::int_vec> true_cuts;
//     CutsGen::rat_vec optLP(3);
//     CutsGen::rat_vec c(3);
//     CutsGen::SpecifiedDirection gd(c);
//     CutsGen::rat_vec b(3);
//
//     A << 1, 0, 0,
//         -5, 1, 0,
//         3, -10, 16;
//
//     c << 1, -5, 3;
//     b << 1, 1, 1;
//
//     Eigen::PartialPivLU<CutsGen::rat_mat> LU(A.transpose().cast<double>());
//
//     optLP = LU.solve(b);
//
//     std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);
//
//     true_cuts.push_back(CutsGen::int_vec(4));
//     true_cuts[0] << 1, 0,  0, std::floor( 143.0/ 16);
//     true_cuts.push_back(CutsGen::int_vec(4));
//     true_cuts[1] << -5, 1, 0, std::floor(-689.0/ 16);
//     true_cuts.push_back(CutsGen::int_vec(4));
//     true_cuts[2] << 3, -10, 1, std::floor(168.0/ 16);
//
//     /*for (const auto& cut : cuts)
//         cout << cut.transpose() << endl;
//
//     for (const auto& cut : true_cuts)
//         cout << cut.transpose() << endl;*/
//
//     assert(cuts == true_cuts);
//     cout << "-----------------------------------------------" << endl;
// }

int main()
{
    try
    {
        seq_min_poly_unit_test();
        pade_approx_unit_test();
        vector_min_poly_unit_test();

        MatrixXi A(3, 3);
        A << 1, 1, 0,
            1, 0, 1,
            0, 1, 1;
        Zp_kernel_unit_test(A, 2);

        A << 2, -1, 3,
            4, 5, 1,
            2, -2, 3;
        Zp_kernel_unit_test(A, 5);
    }
    catch (const std::exception &err)
    {
        ILPCuts::log.out() << err.what() << endl;
        // log.out() << program << endl;
        std::exit(1);
    }
    // catch (const IloException &err)
    // {
    //     log.out() << err << std::endl;
    //     // log.out() << program << endl;
    //     std::exit(1);
    // }
    /*unit_test_I();
    unit_test_S2();
    unit_test_S3();
    unit_test_POW2();*/
}