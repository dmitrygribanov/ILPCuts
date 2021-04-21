// ILPCuts-dev-console.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

#include "conic_cuts_Eigen.h"

//#include "conic_cuts_arageli.h"

//typedef ConicCutsGenerator_Arageli<int, Arageli::rational<int> > CCGen;

int main()
{
    using CutsGen = ILPCuts::ConicCutsGenerator_EigenDence<int, double>;
    CutsGen gen;

    CutsGen::int_mat A, adjA;
    CutsGen::rat_vec c, optLP;
    CutsGen::GreadyDirection gd(c);
    CutsGen::int_type det = 100;

    std::cout

    gen.generate_cuts(A, optLP, gd);
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
