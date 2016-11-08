#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>

template<class Parent>
class dstrat_rwalk_cyclic
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using Parent_t     = typename Parent::This_t;
    using HashSplitter = typename Parent::HashSplitter;
    using Bucket_t     = typename Parent::Bucket_t;

    Parent&      tab;
    std::mt19937 re;
    const size_t steps;
    std::unique_ptr<size_t[]> hist;

    dstrat_rwalk_cyclic(Parent& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps), hist(new size_t[steps])
    {
        //if (! hist) { std::bad_alloc(); }

        for (size_t i = 0; i < steps; ++i)
        {   hist[i] = 0;   }
    }

    dstrat_rwalk_cyclic(Parent& parent, dstrat_rwalk_cyclic&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps), hist(std::move(rhs.hist))
    { }

    bool insert(std::pair<Key,Data> t, HashSplitter hash)
    {
        std::vector<std::pair<std::pair<Key, Data>, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);

        auto tp = t;
        Bucket_t* tb;

        if (bin(re)) tb = tab.getBucket1(hash);
        else         tb = tab.getBucket2(hash);

        queue.emplace_back(tp,tb);
        for (size_t i = 0; !tb->space() && i<steps; ++i)
        {
            auto r = bsd(re);
            tp = tb->replace(r, tp);

            auto hash = tab.h(tp.first);
            if (tab.getBucket1(hash) == tb) tb = tab.getBucket2(hash);
            else                            tb = tab.getBucket1(hash);

            queue.emplace_back(tp,tb);
        }

        if (tb->insert(tp))
        {
            hist[queue.size() -1] += 1;
            return true;
        }

        std::pair<Key,Data> ttp;
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

        return false;
    }
};
