#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template<class Parent>
class dstrat_bfs
{
public:
    using Key            = typename Parent::Key;
    using Data           = typename Parent::Data;
    using Parent_t       = typename Parent::This_t;
    using HashSplitter_t = typename Parent::HashSplitter_t;
    using Bucket_t       = typename Parent::Bucket_t;

    using BFSQueue       = std::vector<std::tuple<Key, int, Bucket_t*> >;

    Parent&      tab;
    const size_t steps;

    dstrat_bfs(Parent& parent, size_t steps = 256, size_t t = 0)
        : tab(parent), steps(steps+1)
    {
        (void)t; /* parameter t is for symmetry with "rwalk" therefore unused*/
    }

    dstrat_bfs(Parent& parent, dstrat_bfs&& rhs)
        : tab(parent), steps(rhs.steps)
    { }

    inline bool expand(BFSQueue& q, size_t index)
    {
        Bucket_t* b = std::get<2>(q[index]);

        for (size_t i = 0; i < tab.bs && q.size() < steps; ++i)
        {
            Key k = b->get(i).first;

            auto hash = tab.h(k);
            Bucket_t* b1 = tab.getBucket1(hash);
            Bucket_t* b2 = tab.getBucket2(hash);

            if        (b == b1)
            {
                q.emplace_back(k, index, b2);
                if (b2->space()) return true;
            }
            else if (b == b2)
            {
                q.emplace_back(k, index, b1);
                if (b1->space()) return true;
            }
            else
            {
                std::cout << "unexpectedly in wrong bucket?" << std::endl;
            }
        }
        return false;
    }

    inline bool rollBackDisplacements(std::pair<Key,Data> t, BFSQueue& bq)
    {
        Key       k1;
        int       prev1;
        Bucket_t* b1;
        std::tie(k1, prev1, b1) = bq[bq.size()-1];

        Key       k2;
        int       prev2;
        Bucket_t* b2;
        while (prev1 >= 0)
        {
            std::tie(k2,prev2,b2) = bq[prev1];

            auto pop = b2->pop(k1);
            if (!pop.first)
            {
                std::cout << "serious issue with rollBack " << k1
                          << " pop " << pop.first << " " << pop.second << std::endl;
                return false;
            }
            if (!b1->insert(k1, pop.second))
            {
                std::cout << "even more serious issue with rollBack" << std::endl;
                return false;
            }

            k1 = k2; prev1 = prev2; b1 = b2;
        }

        if (! b1->insert(t) )
        {   std::cout << "failed final insert" << std::endl; return false; }

        return true;
    }

    inline int insert(std::pair<Key,Data> t, HashSplitter_t hash)
    {
        BFSQueue  bq;
        Bucket_t* b1 = tab.getBucket1(hash);
        Bucket_t* b2 = tab.getBucket2(hash);

        bq.push_back(std::tuple<Key, int, Bucket_t*>(t.first, -1, b1));
        bq.push_back(std::tuple<Key, int, Bucket_t*>(t.first, -1, b2));

        for (size_t i = 0; i < steps; ++i)
        {
            if (expand(bq, i))
            {
                return (rollBackDisplacements(t, bq)) ? bq.size()-2 : -1;
            }
        }

        return -1;
    }
};
