#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_deterministic_walk.hpp
 *
 * dis_deterministic_walk implements an experimental displacement technique
 *similar to random walk, that explores the search space deterministically.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include <utils/debug.hpp>

template <class Parent>
class dis_deterministic_walk_optimistic
{
  private:
    using key_type     = typename Parent::key_type;
    using mapped_type  = typename Parent::mapped_type;
    using value_intern = std::pair<key_type, mapped_type>;

    using parent_type = Parent;
    using hashed_type = typename Parent::hashed_type;
    using bucket_type = typename Parent::bucket_type;

    static constexpr size_t nh = Parent::nh;


    parent_type& tab;
    const size_t steps;

  public:
    dis_deterministic_walk_optimistic(parent_type& parent,
                                      size_t       steps = 256,
                                      size_t             = 0)
        : tab(parent), steps(steps)
    {
    }

    dis_deterministic_walk_optimistic(parent_type& parent,
                                      dis_deterministic_walk_optimistic&& rhs)
        : tab(parent), steps(rhs.steps)
    {
    }

    inline std::pair<int, value_intern*>
    insert(value_intern t, hashed_type hash)
    {
        utils_tm::debug_tm::if_debug(
            "deterministic walk is called on a table with bs>1", tab.bs > 1);

        auto         tp = t;
        bucket_type* tb = tab.get_bucket(hash, 0); // find first bucket
        bucket_type* buckets[nh];

        tp = tb->replace(0, tp); // all buckets are filled thus tp is not dummy
        for (size_t i = 0; i < steps; ++i)
        {
            auto hash = tab.hasher(tp.first);
            tab.get_buckets(hash, buckets); // get all tp's buckets
            size_t which_bucket = 0;        // current bucket of tp + 1
            for (size_t j = 0; j < nh - 1; ++j)
            // this loop just assumes that tp was in nh-1 if it was nowhere else
            {
                if (tb == buckets[j]) which_bucket = j + 1;
            }
            tb = buckets[which_bucket];
            if (tb->space())
            {
                tb->insert(tp);
                return std::make_pair(i, &(tb->elements[0]));
            }
            tp = tb->replace(0, tp);
        }
        return std::make_pair(-1, nullptr);
    }
};
