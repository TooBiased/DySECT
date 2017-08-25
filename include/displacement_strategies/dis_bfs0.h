#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_bfs0.h
 *
 * dis_bfs0 implements the bfs displacement strategy.
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
    class dis_bfs0
    {
    private:
        using key_type       = typename Parent::key_type;
        using mapped_type    = typename Parent::mapped_type;
        using value_intern   = std::pair<key_type,mapped_type>;
        using parent_type    = typename Parent::This_t;
        using hashed_type    = typename Parent::hashed_type;
        using bucket_type    = typename Parent::bucket_type;

        using bfs_queue      = std::vector<std::tuple<key_type, int, bucket_type*> >;

        Parent&      tab;
        const size_t steps;
        static constexpr size_t nh = parent_type::nh;

    public:
        dis_bfs0(Parent& parent, size_t steps = 256, size_t = 0)
            : tab(parent), steps(steps+1)
        {
            /* parameter is for symmetry with "rwalk" therefore unused*/
        }

        dis_bfs0(Parent& parent, dis_bfs0&& rhs)
            : tab(parent), steps(rhs.steps)
        { }

        inline std::pair<int, value_intern*> insert(std::pair<key_type,mapped_type> t, hashed_type hash)
        {
            bfs_queue  bq;

            bucket_type* b[nh];
            tab.get_buckets(hash, b);

            for (size_t i = 0; i < nh; ++i)
            {
                bq.push_back(std::tuple<key_type, int, bucket_type*>(t.first, -1, b[i]));
            }

            for (size_t i = 0; i < steps; ++i)
            {
                if (expand(bq, i))
                {
                    value_intern* pos = rollBackDisplacements(t, bq);
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
                key_type k = b->get(i).first;

                auto hash = tab.hasher(k);

                bucket_type* ptr[nh];
                tab.get_buckets(hash, ptr);
                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (ptr[ti] != b) // POTENTIAL BUG!!! continous bucket problem
                    {
                        //if (overwatch) std::cout << i << " " << ti << std::endl;
                        q.emplace_back(k, index, ptr[ti]);
                        if (ptr[ti]->space()) return true;
                    }
                }
            }
            return false;
        }

        inline value_intern* rollBackDisplacements(std::pair<key_type,mapped_type> t, bfs_queue& bq)
        {
            key_type     k1;
            int          prev1;
            bucket_type* b1;
            std::tie(k1, prev1, b1) = bq[bq.size()-1];

            key_type     k2;
            int          prev2;
            bucket_type* b2;
            while (prev1 >= 0)
            {
                std::tie(k2,prev2,b2) = bq[prev1];

                auto pop = b2->pop(k1);

                b1->insert(k1, pop.second);

                k1 = k2; prev1 = prev2; b1 = b2;
            }

            return b1->insert_ptr(t);
        }
    };

} // namespace cuckoo_displacement
} // namespace dysect
