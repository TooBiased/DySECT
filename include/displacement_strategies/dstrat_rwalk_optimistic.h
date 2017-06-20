#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>

template<class Parent>
class dstrat_rwalk_optimistic
{
private:
    using key_type     = typename Parent::key_type;
    using mapped_type  = typename Parent::mapped_type;
    using value_intern = std::pair<key_type, mapped_type>;

    using Parent_t     = Parent;
    using Hashed_t     = typename Parent::Hashed_t;
    using Bucket_t     = typename Parent::Bucket_t;

    static constexpr size_t nh = Parent::nh;

    Parent_t&    tab;
    std::mt19937 re;
    const size_t steps;

public:
    dstrat_rwalk_optimistic(Parent_t& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps)
    { }

    dstrat_rwalk_optimistic(Parent_t& parent, dstrat_rwalk_optimistic&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
    { }

    inline std::pair<int, value_intern*> insert(value_intern t, Hashed_t hash)
    {
        std::uniform_int_distribution<size_t> bin(0,nh-1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
        //std::uniform_int_distribution<size_t> hfd(0,nh-2);

        auto          tp  = t;
        Bucket_t*     tb  = tab.getBucket(hash, bin(re));
        value_intern* pos = nullptr;

        auto r = bsd(re);
        tp     = tb->replace(r, tp);
        pos    = &(tb->elements[r]);

        for (size_t i = 0; i<steps; ++i)
        {
            auto hash = tab.hasher(tp.first);
            //auto tbd  = tab.getBucket(hash, hfd(re));
            //if (tbd != tb) tb = tbd;
            //else           tb = tab.getBucket(hash, nh-1);
            tb  = tab.getBucket(hash, bin(re));

            if (tb->space()) { tb->insert(tp); return std::make_pair(i, pos); }

            r = bsd(re);
            if (tp.first == t.first) pos = &(tb->elements[r]);
            tp = tb->replace(r, tp);
        }

        return std::make_pair(-1, nullptr);
    }
};
