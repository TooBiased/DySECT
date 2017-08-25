#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_random_walk_optimistic.h
 *
 * dis_random_walk_optimistic implements a random walk displacement
 * technique the displacement path is not stored, therefore, one
 * element will get lost during unsuccessful displacements.
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

template<class Parent>
class dis_random_walk_optimistic
{
private:
    using key_type     = typename Parent::key_type;
    using mapped_type  = typename Parent::mapped_type;
    using value_intern = std::pair<key_type, mapped_type>;

    using parent_type  = Parent;
    using hashed_type  = typename Parent::hashed_type;
    using bucket_type  = typename Parent::bucket_type;

    static constexpr size_t nh = Parent::nh;

    parent_type&    tab;
    std::mt19937 re;
    const size_t steps;

public:
    dis_random_walk_optimistic(parent_type& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re(seed), steps(steps)
    { }

    dis_random_walk_optimistic(parent_type& parent, dis_random_walk_optimistic&& rhs)
        : tab(parent), re(std::move(rhs.re)), steps(rhs.steps)
    { }

    inline std::pair<int, value_intern*> insert(value_intern t, hashed_type hash)
    {
        std::uniform_int_distribution<size_t> bin(0,nh-1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);
        //std::uniform_int_distribution<size_t> hfd(0,nh-2);

        auto          tp  = t;
        bucket_type*  tb  = tab.get_bucket(hash, bin(re));
        value_intern* pos = nullptr;

        auto r = bsd(re);
        tp     = tb->replace(r, tp);
        pos    = &(tb->elements[r]);

        for (size_t i = 0; i<steps; ++i)
        {
            auto hash = tab.hasher(tp.first);
            //auto tbd  = tab.get_bucket(hash, hfd(re));
            //if (tbd != tb) tb = tbd;
            //else           tb = tab.get_bucket(hash, nh-1);
            tb  = tab.get_bucket(hash, bin(re));

            if (tb->space()) { tb->insert(tp); return std::make_pair(i, pos); }

            r = bsd(re);
            if (tp.first == t.first) pos = &(tb->elements[r]);
            tp = tb->replace(r, tp);
        }

        return std::make_pair(-1, nullptr);
    }
};
