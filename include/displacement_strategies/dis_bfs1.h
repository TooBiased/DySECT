#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_bfs1.h
 *
 * dis_bfs1 implements the bfs displacement strategy. This variant
 * is necessary for the overlapping buckets implementation.
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

namespace dysect
{
namespace cuckoo_displacement
{

    template<class Parent>
    class dis_bfs1
    {
    private:
        using key_type       = typename Parent::key_type;
        using mapped_type    = typename Parent::mapped_type;
        using value_intern   = std::pair<key_type,mapped_type>;
        using parent_type    = typename Parent::this_type;
        using hashed_type    = typename Parent::hashed_type;
        using bucket_type    = typename Parent::bucket_type;


        using bfs_item       = std::tuple<value_intern*, int, bucket_type*>;
        using bfs_queue      = std::vector<bfs_item>;

        Parent&      tab;
        const size_t steps;
        static constexpr size_t nh = parent_type::nh;

    public:
        dis_bfs1(Parent& parent, size_t steps = 256, size_t = 0)
            : tab(parent), steps(steps+1)
        { /* parameter is for symmetry with "rwalk" therefore unused*/ }

        dis_bfs1(Parent& parent, dis_bfs1&& rhs)
            : tab(parent), steps(rhs.steps)
        { }

        inline std::pair<int, value_intern*>
        insert(value_intern t, hashed_type hash)
        {
            bfs_queue  bq;
            bucket_type* b[nh];

            tab.get_buckets(hash, b);

            value_intern t_copy = t;

            for (size_t i = 0; i < nh; ++i)
            {
                bq.push_back(bfs_item(&t_copy, -1, b[i]));
            }

            for (size_t i = 0; i < steps; ++i)
            {
                if (expand(bq, i))
                {
                    value_intern* pos = rollBackDisplacements(bq);
                    return std::make_pair((pos) ? bq.size()-nh : -1, pos);
                }
            }

            return std::make_pair(-1, nullptr);
        }

    private:
        inline bool expand(bfs_queue& q, size_t index)
        {
            bucket_type* b = std::get<2>(q[index]);

            for (size_t i = 0; i < tab.bs && q.size() < steps; ++i)
            {
                value_intern* current = &(b->elements[i]);
                key_type      k       = current->first;

                auto hash = tab.hasher(k);

                bucket_type* ptr[nh];
                tab.get_buckets(hash, ptr);
                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (ptr[ti] != b) // POTENTIAL BUG!!! continous bucket problem
                    {
                        //if (overwatch) std::cout << i << " " << ti << std::endl;
                        q.emplace_back(current, index, ptr[ti]);
                        if (ptr[ti]->space()) return true;
                    }
                }
            }
            return false;
        }

        inline value_intern* rollBackDisplacements(bfs_queue& bq)
        {
            value_intern* k1;
            int           prev1;
            bucket_type*  b1;
            std::tie(k1, prev1, b1) = bq[bq.size()-1];

            value_intern* t = b1->probe_ptr(key_type()).second;
            (*t) = *k1;

            value_intern* k2;
            int           prev2;
            bucket_type*  b2;

            value_intern* k0 = nullptr;
            while (prev1 >= 0)
            {
                std::tie(k2,prev2,b2) = bq[prev1];
                (*k1) = *k2;

                k0 = k1;
                k1 = k2; prev1 = prev2; b1 = b2;
            }

            return k0;
        }

    };

} // namespace cuckoo_displacement
} // namespace dysect
