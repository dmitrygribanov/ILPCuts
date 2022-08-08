#pragma once

#include <fstream>
#include <vector>

namespace ILPCuts
{
    namespace Utility
    {
        namespace primes
        {
            class DownloadException : public Exception
            {
            public:
                DownloadException()
                    : Exception("ILPCuts error: could not download primes table...") {}

                DownloadException(string path)
                    : Exception("ILPCuts error: could not download primes table; path = \"" + path + "\"") {}
            };

            const std::string default_path = "../primes_table/primes_mill.txt";
            const size_t max_download = 1000000;

            template <typename ptype>
            std::vector<ptype> download(size_t max_num, std::string path = default_path)
            {
                std::ifstream file(path);

                if (!file)
                    throw DownloadException(path);

                std::vector<ptype> res;

                size_t i = 0;
                while (file && (i++ < max_num))
                {
                    res.push_back(0);
                    file >> res.back();
                }

                return std::move(res);
            }

            template <typename int_type, typename lint_type>
            int_type find_divisor(const lint_type &val, const std::vector<int_type> &primes_data)
            {
                for (int_type p : primes_data)
                {
                    if (val % p == 0)
                        return p;
                }

                return 0;
            }
        }
    } // namespace Utility

} // namespace ILPCuts
