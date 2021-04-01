// Tests.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

#include "conic_cuts.h"

#include <iostream>

#include <arageli.hpp>

void unit_test_I() {
    typedef ConicCutsGenerator_Arageli<int, Arageli::rational<int> > CCGen;

    CCGen cc_gen;
    CCGen::int_mat Ab = "((1,0,0,0,1),\
        (0,1,0,0,2),\
        (1,1,0,0,3),\
        (1,-1,0,0,-1),\
        (0,0,1,0,4),\
        (0,0,0,1,-5))";
    CCGen::rational_vec vertex = "(1,2,4,-5)";
    CCGen::int_vec cT = "(1,1,1,0)";
    CCGen::int_mat cuts = cc_gen.make_cuts(vertex, cT, Ab);

    assert(cuts == CCGen::int_mat("((1,0,0,0,1),\
        (0,1,0,0,2),\
        (0,0,1,0,4),\
        (0,0,0,1,-5))")
    );
}

void unit_test_S2() {
    typedef ConicCutsGenerator_Arageli<int, Arageli::rational<int> > CCGen;

    CCGen cc_gen;
    CCGen::int_mat Ab = "((1,0,0,1),(0,1,0,2),(0,0,2,11),(1,-1,0,-1))";
    CCGen::rational_vec vertex = "(1, 2, 11/2)";
    CCGen::int_vec cT = "(1,1,1)";
    CCGen::int_mat cuts = cc_gen.make_cuts(vertex, cT, Ab);

    assert(cuts == CCGen::int_mat("((1,0,0,1),\
        (0,1,0,2),\
        (0,0,1,5)")
    );
}

void unit_test_S3() {
    typedef ConicCutsGenerator_Arageli<int, Arageli::rational<int> > CCGen;

    CCGen cc_gen;
    CCGen::int_mat Ab = "((1,0,0,1),(0,1,0,2),(0,0,3,11),(1,-1,0,-1))";
    CCGen::rational_vec vertex = "(1, 2, 11/3)";
    CCGen::int_vec cT = "(1,1,1)";
    CCGen::int_mat cuts = cc_gen.make_cuts(vertex, cT, Ab);

    assert(cuts == CCGen::int_mat("((1,0,0,1),\
        (0,1,0,2),\
        (0,0,1,3)")
    );
}

int main()
{
    unit_test_I();
    unit_test_S2();
    unit_test_S3();

    std::cout << "Hello World!\n";
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
