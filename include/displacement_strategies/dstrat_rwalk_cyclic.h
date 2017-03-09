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
    using Pair_t       = std::pair<Key, Data>;
    using Parent_t     = Parent;
    using Hashed_t     = typename Parent::Hashed_t;
    using Bucket_t     = typename Parent::Bucket_t;

    static constexpr size_t nh = Parent::nh;

    Parent_t&    tab;
    std::mt19937 re;
    const size_t steps;

    dstrat_rwalk_cyclic(Parent_t& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps)
    { }

    dstrat_rwalk_cyclic(Parent_t& parent, dstrat_rwalk_cyclic&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
    { }

    inline std::pair<int, Pair_t*> insert(Pair_t t, Hashed_t hash)
    {
        std::vector<std::pair<Pair_t, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,nh-1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
        std::uniform_int_distribution<size_t> hfd(0,nh-2);

        auto tp = t;
        Bucket_t* tb = tab.getBucket(hash, bin(re));
        Pair_t* pos = nullptr;

        queue.emplace_back(tp,tb);
        for (size_t i = 0; !tb->space() && i<steps; ++i)
        {
            auto r = bsd(re);
            if (tp.first == t.first) pos = &(tb->elements[r]);
            tp = tb->replace(r, tp);

            auto hash = tab.hasher(tp.first);
            auto tbd  = tab.getBucket(hash, hfd(re));
            if (tbd != tb) tb = tbd;
            else           tb = tab.getBucket(hash, nh-1);

            queue.emplace_back(tp,tb);
        }

        if (tb->insert(tp))
        {
            return std::make_pair(queue.size() -1, pos);
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

        return std::make_pair(-1, nullptr);
    }
};
