#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_walk_acyclic.h
 *
 * dis_random_walk_acyclic
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

    template<class Parent> // think about spezialization if nh = 2
    class dis_random_walk_acyclic
    {
    private:
        using key_type       = typename Parent::key_type;
        using mapped_type    = typename Parent::mapped_type;
        using value_intern   = std::pair<key_type,mapped_type>;

        using parent_type    = Parent;
        using hashed_type    = typename Parent::hashed_type;
        using bucket_type    = typename Parent::bucket_type;

        parent_type& tab;
        std::mt19937 re;
        const size_t steps;

        static constexpr size_t nh = Parent::nh;

    public:
        dis_random_walk_acyclic(parent_type& parent, size_t steps=256, size_t seed=30982391937209388ull)
            : tab(parent), re(seed), steps(steps)
        { }

        dis_random_walk_acyclic(parent_type& parent, dis_random_walk_acyclic&& rhs)
            : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
        { }

        inline std::pair<int, value_intern*> insert(std::pair<key_type,mapped_type> t, hashed_type hash)
        {

            std::vector<std::pair<value_intern, bucket_type*> > queue;
            std::uniform_int_distribution<size_t> bin(0,nh-1);
            std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
            std::uniform_int_distribution<size_t> hfd(0,nh-2);

            auto         tp = t;
            bucket_type* tb = tab.get_bucket(hash, bin(re));

            queue.emplace_back(tp,tb);

            size_t i = 0;
            for ( ; !(tb->space()) && i<steps; ++i)
            {
                auto r = bsd(re);
                tp = tb->get(r);
                auto hash = tab.hasher(tp.first);
                auto tbd = tab.get_bucket(hash,hfd(re));
                if (tbd != tb) tb = tbd;
                else           tb = tab.get_bucket(hash, nh-1);

                queue.emplace_back(tp,tb);

                // Explicit Cycle Detection (they are popped from the displacement queue)
                for (size_t j = 0; j < queue.size() - 1; ++j)
                {
                    if (queue[j].second == tb)
                    {
                        while (queue.size() > j+1) queue.pop_back();
                        break;
                    }
                }
            }

            if (! tb->space()) { return std::make_pair(-1, nullptr); }

            for (size_t i = queue.size()-1; i > 0; --i)
            {
                std::tie(tp,tb) = queue[i];
                if (! queue[i-1].second->remove(tp.first) ||
                    ! tb->insert(tp.first, tp.second))
                { return std::make_pair(-1, nullptr); }

            }

            value_intern* pos = queue[0].second->insert_ptr(t);

            return std::make_pair((pos) ? i : -1, pos);
        }
    };
} // namespace cuckoo_displacement
} // namespace dysect
