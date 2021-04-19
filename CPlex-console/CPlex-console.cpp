
#include <iostream>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "conic_cuts_Eigen.h"

ILOUSERCUTCALLBACK0(myCutsCallback) {
    IloEnv env = getEnv();
    IloModel mod = getModel();
    env.out() << "--------------My cut Callback----------------" << std::endl;
    IloModel::Iterator iter = IloModel::Iterator(mod);
    while (iter.ok()) {
        env.out() << (*iter).asNumExpr() << std::endl;
        ++iter;
    }
    env.out() << "---------------------------------------------" << std::endl;
}

ILOMIPINFOCALLBACK0(loggingCallback)
{
   IloEnv env = getEnv();
   IloModel mod = getModel();
   env.out() << "--------------My info Callback---------------" << std::endl;
   env.out() << "---------------------------------------------" << std::endl;
}

int main(int argc, char** argv)
{
    IloEnv env;
    try {
        IloModel model(env);
        IloCplex cplex(env);
        IloObjective   obj;
        IloNumVarArray vars(env);
        IloRangeArray  rngs(env);

        const char* datadir = (argc >= 2) ? argv[1] : "../problems/cplex";
        stringstream filename;
        filename << datadir << "/" << "noswot.mps";
        env.out() << "reading " << filename.str() << endl;
        cplex.importModel(model, filename.str().c_str(), obj, vars, rngs);

        env.out() << "extracting model ..." << endl;
        cplex.extract(model);

        cplex.setParam(IloCplex::Param::MIP::Interval, 1000);
        env.out() << "solving model ...\n";

        cplex.use(myCutsCallback(env));

        if (!cplex.solve()) {
            env.error() << "Failed to optimize." << endl;
            throw(-1);
        }

        env.out() << "solution status is " << cplex.getStatus() << endl;
        env.out() << "solution value  is " << cplex.getObjValue() << endl;
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }

    env.end();

    return 0;
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
