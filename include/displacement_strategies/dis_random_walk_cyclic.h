#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_random_walk_cyclic.h
 *
 * dis_random_walk_cyclic implements a random walk displacement
 * technique, which precomputes the whole displacement path. It does
 * not remove possible cycles, instead cycles get handled during the
 * actual displacement.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <tuple>
#include <random>

namespace dysect
{
namespace cuckoo_displacement
{

    template<class Parent>
    class dis_random_walk_cyclic
    {
    private:
        using key_type     = typename Parent::key_type;
        using mapped_type  = typename Parent::mapped_type;
        using value_intern = std::pair<key_type, mapped_type>;

        using parent_type  = Parent;
        using hashed_type  = typename Parent::hashed_type;
        using bucket_type  = typename Parent::bucket_type;

        static constexpr size_t nh = Parent::nh;

        parent_type& tab;
        std::mt19937 re;
        const size_t steps;

    public:
        dis_random_walk_cyclic(parent_type& parent, size_t steps=256, size_t seed=30982391937209388ull)
            : tab(parent), re(seed), steps(steps)
        { }

        dis_random_walk_cyclic(parent_type& parent, dis_random_walk_cyclic&& rhs)
            : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
        { }

        inline std::pair<int, value_intern*> insert(value_intern t, hashed_type hash)
        {
            std::vector<std::pair<value_intern, bucket_type*> > queue;
            std::uniform_int_distribution<size_t> bin(0,nh-1);
            std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
            std::uniform_int_distribution<size_t> hfd(0,nh-2);

            auto tp = t;
            bucket_type*  tb  = tab.get_bucket(hash, bin(re));
            value_intern* pos = nullptr;

            queue.emplace_back(tp,tb);
            for (size_t i = 0; !tb->space() && i<steps; ++i)
            {
                auto r = bsd(re);
                if (tp.first == t.first) pos = &(tb->elements[r]);
                tp = tb->replace(r, tp);

                auto hash = tab.hasher(tp.first);
                auto tbd  = tab.get_bucket(hash, hfd(re));
                if (tbd != tb) tb = tbd;
                else           tb = tab.get_bucket(hash, nh-1);

                queue.emplace_back(tp,tb);
            }

            if (tb->insert(tp))
            {
                return std::make_pair(queue.size() -1, pos);
            }

            std::pair<key_type,mapped_type> ttp;
            std::tie(ttp, tb) = queue[queue.size() - 1];
            for (size_t i = queue.size() - 2; i >= 1; --i)
            {
                std::tie(tp, tb) = queue[i];
                if (!(tb->remove(tp.first)))
                { std::cout << "f1" << std::endl; }
                if (!(tb->insert(ttp.first, ttp.second)))
                { std::cout << "f2" << std::endl; };
                ttp = tp;
            }

            return std::make_pair(-1, nullptr);
        }
    };

} // namespace cuckoo_displacement
} // namespace dysect
