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

template <class Parent>
class dstrat_triv
{
public:
    using key_type            = typename Parent::key_type;
    using mapped_type         = typename Parent::mapped_type;
    using Hashed_t            = typename Parent::Hashed_t;

    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    inline int insert(std::pair<key_type, mapped_type>, Hashed_t)
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
struct Config
{
    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;
    static constexpr size_t nh = NH;
    static constexpr size_t sbs = 4;

    template <class T>
    using DisStrat_temp = DisStrat<T>;

    using HistCount_t = HistCount;
};

template<size_t NS = 64, typename GRat = std::ratio<11,10> >
struct HopscotchConfig
{
    static constexpr size_t NeighborSize = NS;
    static constexpr double GrowRatio_d  = double(GRat::num)/double(GRat::den);
    using GrowRatio = GRat;
};

template<class HistCount = no_hist_count>
struct TrivConfig
{
    using HistCount_t = HistCount;
};
