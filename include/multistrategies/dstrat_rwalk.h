#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>



template<class Parent> // think about spezialization if nh = 2
class dstrat_multi_rwalk
{
public:
    using Key            = typename Parent::Key;
    using Data           = typename Parent::Data;
    using Parent_t       = Parent;
    using HashSplitter_t = typename Parent::HashSplitter_t;
    using Bucket_t       = typename Parent::Bucket_t;

    Parent&      tab;
    std::mt19937 re;
    const size_t steps;

    static constexpr size_t nh = Parent::nh;

    dstrat_multi_rwalk(Parent_t& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps)
    { }

    dstrat_multi_rwalk(Parent_t& parent, dstrat_multi_rwalk&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
    { }

    inline int insert(std::pair<Key,Data> t, HashSplitter_t hash)
    {
        std::vector<std::pair<std::pair<Key, Data>, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,nh-1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
        std::uniform_int_distribution<size_t> hfd(0,nh-2);

        auto      tp = t;
        Bucket_t* tb = tab.getBucket(hash, bin(re));

        queue.emplace_back(tp,tb);

        size_t i = 0;
        for ( ; !(tb->space()) && i<steps; ++i)
        {
            auto r = bsd(re);
            tp = tb->get(r);
            auto hash = tab.h(tp.first);
            auto tbd = tab.getBucket(hash,hfd(re));
            if (tbd != tb) tb = tbd;
            else           tb = tab.getBucket(hash, nh-1);

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

        if (! tb->space()) { return -1; }

        //hist[queue.size() - 1] += 1;

        for (size_t i = queue.size()-1; i > 0; --i)
        {
            std::tie(tp,tb) = queue[i];
            if (! queue[i-1].second->remove(tp.first))  { std::cout << "e2" << std::endl; return -1; }
            if (! tb->insert(tp.first, tp.second))      { std::cout << "e1" << std::endl; return -1; }
        }

        if (! queue[0].second->insert(t)) { std::cout << "e3" << std::endl; return -1; }

        return i;
    }
};
