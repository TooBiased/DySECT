#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template<class Parent>
class dstrat_bfs
{
private:
    using key_type       = typename Parent::key_type;
    using mapped_type    = typename Parent::mapped_type;
    using value_intern   = std::pair<key_type,mapped_type>;
    using Parent_t       = typename Parent::This_t;
    using Hashed_t       = typename Parent::Hashed_t;
    using Bucket_t       = typename Parent::Bucket_t;

    using BFSQueue       = std::vector<std::tuple<key_type, int, Bucket_t*> >;

    Parent&      tab;
    const size_t steps;
    static constexpr size_t nh = Parent_t::nh;

public:
    dstrat_bfs(Parent& parent, size_t steps = 256, size_t = 0)
        : tab(parent), steps(steps+1)
    {
        /* parameter is for symmetry with "rwalk" therefore unused*/
    }

    dstrat_bfs(Parent& parent, dstrat_bfs&& rhs)
        : tab(parent), steps(rhs.steps)
    { }

    inline std::pair<int, value_intern*> insert(std::pair<key_type,mapped_type> t, Hashed_t hash)
    {
        BFSQueue  bq;

        Bucket_t* b[nh];
        tab.getBuckets(hash, b);

        for (size_t i = 0; i < nh; ++i)
        {
            bq.push_back(std::tuple<key_type, int, Bucket_t*>(t.first, -1, b[i]));
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
    inline bool expand(BFSQueue& q, size_t index)
    {
        Bucket_t* b = std::get<2>(q[index]);

        for (size_t i = 0; i < tab.bs && q.size() < steps; ++i)
        {
            key_type k = b->get(i).first;

            auto hash = tab.hasher(k);

            Bucket_t* ptr[nh];
            tab.getBuckets(hash, ptr);
            for (size_t ti = 0; ti < nh; ++ti)
            {
                if (ptr[ti] != b)
                {
                    //if (overwatch) std::cout << i << " " << ti << std::endl;
                    q.emplace_back(k, index, ptr[ti]);
                    if (ptr[ti]->space()) return true;
                }
            }
        }
        return false;
    }

    inline value_intern* rollBackDisplacements(std::pair<key_type,mapped_type> t, BFSQueue& bq)
    {
        key_type       k1;
        int       prev1;
        Bucket_t* b1;
        std::tie(k1, prev1, b1) = bq[bq.size()-1];

        key_type       k2;
        int       prev2;
        Bucket_t* b2;
        while (prev1 >= 0)
        {
            std::tie(k2,prev2,b2) = bq[prev1];

            auto pop = b2->pop(k1);

            b1->insert(k1, pop.second);

            k1 = k2; prev1 = prev2; b1 = b2;
        }

        return b1->insertPtr(t);
    }
};
