// Tests.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

#include "conic_cuts_Eigen.h"

//#include "conic_cuts_arageli.h"
using std::cout;
using std::endl;

void unit_test_I() {

    using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<int, double>;
    CutsGen gen;

    CutsGen::int_mat        A(4,4);
    CutsGen::rat_vec      optLP(4);
    CutsGen::rat_vec          c(4);
    CutsGen::GreadyDirection gd(c);

    A << 1,   0,   0, 0,
         1,   1,   0, 0,
        10, -10,   1, 0,
        24,  30, -30, 1;

    CutsGen::int_vec t(4);
    t << 1, 1, 1, 1;

    c = (A * t).cast<double>();

    optLP << 1, 2, 4, -5;

    std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);

    cout << "---------------unit_test_I---------------------" << endl;
    assert(cuts.empty());

    cout << gen.out().str() << endl;
    cout << gen.err().str() << endl;
    cout << "-----------------------------------------------" << endl;
}

void unit_test_S2() {
    using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<int, double>;
    CutsGen gen;

    CutsGen::int_mat A(3, 3);
    std::vector<CutsGen::int_vec> true_cuts;
    CutsGen::rat_vec optLP(3);
    CutsGen::rat_vec c(3);
    CutsGen::GreadyDirection gd(c);

    A << 1, 0, 0,
         1, 1, 0,
         0, 0, 2;

    c << 1, 1, 1;

    optLP << 1, 2, 5.5;

    std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);
    
    true_cuts.push_back(CutsGen::int_vec(4));

    true_cuts[0] << 0, 0, 1, 5;

    cout << "---------------unit_test_S2---------------------" << endl;
    assert( cuts == true_cuts );

    cout << gen.out().str() << endl;
    cout << gen.err().str() << endl;
    cout << "-----------------------------------------------" << endl;
}

void unit_test_S3() {
    using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<int, double>;
    CutsGen gen;

    CutsGen::int_mat A(3, 3);
    std::vector<CutsGen::int_vec> true_cuts;
    CutsGen::rat_vec optLP(3);
    CutsGen::rat_vec c(3);
    CutsGen::GreadyDirection gd(c);

    A << 1, 0, 0,
        1, 1, 0,
        0, 0, 3;

    c << 1, 1, 1;

    optLP << 1, 2, (float)(11) / 3;

    std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);

    true_cuts.push_back(CutsGen::int_vec(4));

    true_cuts[0] << 0, 0, 1, 3;

    cout << "---------------unit_test_S3---------------------" << endl;
    assert(cuts == true_cuts);

    cout << gen.out().str() << endl;
    cout << gen.err().str() << endl;
    cout << "-----------------------------------------------" << endl;
}

void unit_test_POW2() {
    using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<int, double>;
    CutsGen gen;

    CutsGen::int_mat A(3, 3);
    std::vector<CutsGen::int_vec> true_cuts;
    CutsGen::rat_vec optLP(3);
    CutsGen::rat_vec c(3);
    CutsGen::GreadyDirection gd(c);
    CutsGen::rat_vec b(3);

    A << 1, 0, 0,
        -5, 1, 0,
        3, -10, 16;

    c << 1, -5, 3;
    b << 1, 1, 1;

    Eigen::PartialPivLU<CutsGen::rat_mat> LU(A.transpose().cast<double>());

    optLP = LU.solve(b);

    std::vector<CutsGen::int_vec> cuts = gen.generate_cuts(A, optLP, gd);

    true_cuts.push_back(CutsGen::int_vec(4));
    true_cuts[0] << 1, 0,  0, std::floor( 143.0/ 16);
    true_cuts.push_back(CutsGen::int_vec(4));
    true_cuts[1] << -5, 1, 0, std::floor(-689.0/ 16);
    true_cuts.push_back(CutsGen::int_vec(4));
    true_cuts[2] << 3, -10, 1, std::floor(168.0/ 16);

    cout << "---------------unit_test_POW2---------------------" << endl;

    /*for (const auto& cut : cuts)
        cout << cut.transpose() << endl;

    for (const auto& cut : true_cuts)
        cout << cut.transpose() << endl;*/

    cout << gen.out().str() << endl;
    cout << gen.err().str() << endl;

    assert(cuts == true_cuts);
    cout << "-----------------------------------------------" << endl;
}

int main()
{
    unit_test_I();
    unit_test_S2();
    unit_test_S3();
    unit_test_POW2();
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
