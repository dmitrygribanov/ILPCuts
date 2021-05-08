// ILPCuts-dev-console.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <mpir.h>
#include <mpirxx.h>

#include "mpir_eigen_supp.h"
#include "conic_cuts_Eigen.h"


int main()
{
    using CutsGen = ILPCuts::ConicCutsGenerator<mpz_class, mpq_class>;
    CutsGen gen;
    gen.Init();

    CutsGen::is_mat A, adjA;
    CutsGen::rs_vec c, optLP;
    //CutsGen::GreadyDirection gd(c);
    CutsGen::int_type det = 100;

    gen.generate_cuts(A, optLP);

    ///*mpf_class x("1234444444441354231341432415234.0000099999999999");
    //mpf_class y("12344444434534534526246.9999999999999");*/

    //mpf_class x("1234444444441354231341432415234.0000099999999999", 10000);
    //mpf_class y("123444444444444444444444444444444444444444444.999999999999999999999999999999999", 10000);

    //mpz_class xi = static_cast<mpz_class>(round(x));
    //mpz_class yi = static_cast<mpz_class>(round(y));

    //std::cout << x << ", " << xi << ", " << abs(x-xi) << std::endl;
    //std::cout << y << ", " << yi << ", " << abs(y - yi) << std::endl;

	/*if (x > mpf_class(0)) {
		y = floor(x + mpf_class(0.5));
	}
	else {
		y = ceil(x - mpf_class(0.5));
	}*/
    
    //y = (x > mpf_class(0)) ? floor(x + mpf_class(0.5)) : ceil(x - mpf_class(0.5));
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
