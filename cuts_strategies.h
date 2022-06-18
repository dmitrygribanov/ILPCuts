#pragma once

#include "basic_utils.h"

namespace ILPCuts
{
    template <typename intT, typename ratT>
    struct CutsGenStrategies
    {
        using int_type = intT;
        using rational_type = ratT;
        using is_mat = SparseMatrix<int_type>;
        using is_vec = SparseVector<int_type>;
        using rs_vec = SparseVector<rational_type>;

        class DirectionsStrategyBase
        {
        public:
            virtual std::vector<rs_vec> generate_directons(const is_mat &A) = 0;
        };

        class GreadyDirection : public DirectionsStrategyBase
        {
        private:
            rs_vec _c;

        public:
            GreadyDirection(const rs_vec &c)
                : _c(c) {}

            std::vector<rs_vec> generate_directons(const is_mat &A) override
            {
                std::vector<rs_vec> res(1);
                res[0] = _c;
                return std::move(res);
            }
        };

        class SomeDirection : public DirectionsStrategyBase
        {
        public:
            std::vector<rs_vec> generate_directons(const is_mat &A) override
            {
                size_t n = A.cols();
                rs_vec direction(n);
                direction.setZero();

                for (size_t i = 0; i < n; ++i)
                {
                    direction += A.template cast<rational_type>().col(i);
                }
                direction /= rational_type(n);

                std::vector<rs_vec> res(1);
                res[0] = std::move(direction);
                return std::move(res);
            }
        };

        template <typename randGenT>
        class ConvUniform : public DirectionsStrategyBase
        {
        public:
            using gen_type = randGenT;

        private:
            const size_t _rand_dir_num;
            gen_type &_gen;

        public:
            ConvUniform(gen_type &gen, size_t rand_dir_num)
                : DirectionsStrategyBase(), _gen(gen), _rand_dir_num(rand_dir_num) {}

            std::vector<rs_vec> generate_directons(const is_mat &A) override
            {
                std::vector<rs_vec> result_directions;
                size_t n = A.cols();

                result_directions.reserve(_rand_dir_num);

                std::uniform_real_distribution<double> udistr(0.0, 1.0);

                for (size_t i = 0; i < _rand_dir_num; ++i)
                {
                    rs_vec new_dir = A.col(0).template cast<rational_type>();
                    for (size_t j = 1; j < n; ++j)
                    {
                        rational_type alpha = udistr(_gen);
                        new_dir = (new_dir * alpha + A.col(j).template cast<rational_type>() * (1.0 - alpha)).eval();
                    }
                    result_directions.push_back(std::move(new_dir));
                }

                return std::move(result_directions);
            }
        };
    };

}