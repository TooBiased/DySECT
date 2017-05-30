#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>



template<class Parent> // think about spezialization if nh = 2
class dstrat_rwalk_anticycle
{
private:
    using key_type       = typename Parent::key_type;
    using mapped_type    = typename Parent::mapped_type;
    using value_intern   = std::pair<key_type,mapped_type>;

    using Parent_t       = Parent;
    using Hashed_t       = typename Parent::Hashed_t;
    using Bucket_t       = typename Parent::Bucket_t;

    Parent_t&    tab;
    std::mt19937 re;
    const size_t steps;

    static constexpr size_t nh = Parent::nh;

public:
    dstrat_rwalk_anticycle(Parent_t& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps)
    { }

    dstrat_rwalk_anticycle(Parent_t& parent, dstrat_rwalk_anticycle&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
    { }

    inline std::pair<int, value_intern*> insert(std::pair<key_type,mapped_type> t, Hashed_t hash)
    {

        std::vector<std::pair<value_intern, Bucket_t*> > queue;
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
            auto hash = tab.hasher(tp.first);
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

        if (! tb->space()) { return std::make_pair(-1, nullptr); }

        for (size_t i = queue.size()-1; i > 0; --i)
        {
            std::tie(tp,tb) = queue[i];
            if (! queue[i-1].second->remove(tp.first) ||
                ! tb->insert(tp.first, tp.second))
            { return std::make_pair(-1, nullptr); }

        }

        value_intern* pos = queue[0].second->insertPtr(t);

        return std::make_pair((pos) ? i : -1, pos);
    }
};
