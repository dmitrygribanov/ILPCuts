#include "wiedemann_alg.h"

using namespace ILPCuts;
using namespace std;

namespace ILPCuts
{
    logger log("tests_log.txt");

    void seq_min_poly_unit_test()
    {
        log.out() << "Sequence minimal polynomial test..." << endl;

        zz_p::init(5);

        stringstream str;
        zz_pX fX;
        str << "[3 0 4 2 3 0]";
        str >> fX;

        log.out() << "sequence: " << str.str() << endl;

        zz_pX mX = seq_min_poly(fX);
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

        stringstream str;
        zz_pX fX;
        str << "[1 2 3 4]";
        str >> fX;

        log.out() << "sequence(x) = " << str.str() << "\n"
                  << "(k, n-k) = (2,2)" << endl;

        zz_pX sX, tX;
        tie(sX, tX) = pade_approx(fX, 2);

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
//     CutsGen::GreadyDirection gd(c);
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
//     CutsGen::GreadyDirection gd(c);
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
//     CutsGen::GreadyDirection gd(c);
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
//     CutsGen::GreadyDirection gd(c);
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