#pragma once

/*******************************************************************************
 * include/config.h
 *
 * collection of configuration types that are mostly used to reduce
 * the number of necessary template parameters, used by our hash
 * tables.  This reduces the exposure of implementation details to the
 * user.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <memory>
#include <iostream>
#include <tuple>
#include <ratio>

namespace dysect
{

    template <class Parent>
    class dstrat_triv
    {
    public:
        using key_type            = typename Parent::key_type;
        using mapped_type         = typename Parent::mapped_type;
        using hashed_type         = typename Parent::hashed_type;
        using value_intern        = std::pair<key_type, mapped_type>;

        dstrat_triv(Parent&, size_t, size_t) {}
        dstrat_triv(Parent&, dstrat_triv&&) {}

        inline int insert(value_intern, hashed_type)
        {   return -1; }
    };


    class hist_count
    {
    public:
        hist_count(size_t s) : steps(s), hist(new size_t[s])
        { for (size_t i = 0; i < s; ++i) { hist[i] = 0; } }

        void add(size_t i) { auto ind = (i<steps) ? i:steps-1; ++hist[ind];}

        const size_t steps;
        std::unique_ptr<size_t[]> hist;
    };


    class no_hist_count
    {
    public:
        no_hist_count(size_t = 0) { }
        void add(size_t) { }
        static constexpr size_t  steps = 0;
        static constexpr size_t* hist  = nullptr;
    };

    template<size_t BS = 8, size_t NH = 2, size_t TL = 256,
             template <class> class DisStrat = dstrat_triv,
             class HistCount = hist_count>
    struct cuckoo_config
    {
        static constexpr size_t bs = BS;
        static constexpr size_t tl = TL;
        static constexpr size_t nh = NH;
        static constexpr size_t sbs = 4;

        template <class T>
        using dis_strat_type  = DisStrat<T>;

        using hist_count_type = HistCount;
    };

    template<size_t NS = 64, typename GRat = std::ratio<11,10> >
    struct hopscotch_config
    {
        static constexpr size_t neighborhood_size = NS;
        static constexpr double grow_ratio_double = double(GRat::num)/double(GRat::den);
        using grow_ratio = GRat;
    };

    template<class HistCount = no_hist_count>
    struct triv_config
    {
        using hist_count_type = HistCount;
    };

} // namespace dysect
