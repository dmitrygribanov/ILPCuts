#pragma once

#include "cuts_gen_cone_coverer.h"

namespace ILPCuts
{
    class ConstantSizeWindow
    {
    public:
        size_t wSize;

        ConstantSizeWindow(size_t inputSize)
            : wSize(inputSize) {}

        size_t operator()(const local_instance &) const { return wSize; }
    };

    class DimSizeWindow
    {
    public:
        size_t maxSize = 0;

        DimSizeWindow(size_t aMaxSize)
            : maxSize(aMaxSize)
        {
        }

        size_t operator()(const local_instance &inst) const
        {
            size_t n = inst.get_num_vars();

            if (n > maxSize)
                return maxSize;

            return n;
        }
    };

    class SqrtDimSizeWindow
    {
    public:
        size_t maxSize = 0;
        size_t minSize = 0;

        SqrtDimSizeWindow(size_t aMinSize, size_t aMaxSize)
            : maxSize(aMaxSize), minSize(aMinSize)
        {
        }

        size_t operator()(const local_instance &inst) const
        {
            size_t n = inst.get_num_vars();
            double sqrt_n = std::floor(
                std::sqrt(static_cast<double>(n)));

            if (sqrt_n < minSize)
                return minSize;

            if (sqrt_n > maxSize)
                return maxSize;

            return sqrt_n;
        }
    };

    template <typename rgenT>
    class RandomNondegSubmatGen
    {
    public:
        size_t numberOfTrials = 5;
        std::mt19937 *rgen = nullptr;

        RandomNondegSubmatGen() {}

        template <typename submatSizeFunT>
        std::tuple<std::vector<size_t>, std::vector<size_t>> operator()(const local_instance &inst, submatSizeFunT submat_size_fun) const
        {
            size_t n = inst.get_num_vars();
            size_t submat_n = submat_size_fun(n);

            std::vector<size_t> idxs(n);
            std::iota(idxs.begin(), idxs.end(), 0);

            for (size_t trial = 0; trial < numberOfTrials; ++trial)
            {
                std::vector cols, rows;
                std::sample(idxs.begin(), idxs.end(), std::back_inserter(rows), submat_n, *rgen);
                std::sample(idxs.begin(), idxs.end(), std::back_inserter(cols), submat_n, *rgen);
                std::sort(rows.begin(), rows.end());
                std::sort(cols.begin(), cols.end());

                auto B = generate_square_submatrix_byidx(inst.A, rows, cols);

                if (rank(B) = submat_n)
                    return std::make_tuple(rows, cols);
            }
        }
    };

} // namespace ILPCuts
