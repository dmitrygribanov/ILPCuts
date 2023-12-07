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

        double end() const
        {
            return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - _start).count();
        }

        string finish_str(string suff = "'s\n") const
        {
            stringstream res;

            res << end() << suff;

            return res.str();
        }
    };

    Timer global_timer;
} // namespace ILPCuts
