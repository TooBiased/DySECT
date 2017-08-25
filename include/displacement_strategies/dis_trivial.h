#pragma once

/*******************************************************************************
 * include/displacement_strategies/dis_trivial.h
 *
 * dstrat_trivial is the trivialial displacement technique (no displacements
 * at all). The hash table will be similar to k-choice bucket hashing.
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

    template <class Parent>
    class dis_trivial
    {
    private:
        using key_type      = typename Parent::key_type;
        using mapped_type   = typename Parent::mapped_type;
        using value_intern  = std::pair<key_type,mapped_type>;

        using hashed_type   = typename Parent::hashed_type;

    public:
        dis_trivial(Parent&, size_t, size_t) {}
        dis_trivial(Parent&, dis_trivial&&) {}

        inline std::pair<int, value_intern*> insert(value_intern, hashed_type)
        {   return std::make_pair(-1, nullptr); }
    };

} // namespace cuckoo_displacement
} // namespace dysect
