#pragma once
#include <chrono>
#include <sstream>

namespace ILPCuts
{
    using namespace std;

    class Timer
    {
    private:
        mutable chrono::steady_clock::time_point _start;

    public:
        void start() const
        {
            _start = chrono::steady_clock::now();
        }

        string finish_str(string suff = "'s\n") const
        {
            stringstream res;

            auto end = std::chrono::steady_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - _start;

            res << elapsed_seconds.count() << suff;

            return res.str();
        }
    };

    Timer global_timer;
} // namespace ILPCuts
