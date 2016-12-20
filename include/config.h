#pragma once

#include <memory>
#include "include/multistrategies/dstrat_triv.h"
#include "include/multistrategies/dstrat_bfs.h"

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
    no_hist_count(size_t) { }
    void add(size_t) { }
    static constexpr size_t  steps = 0;
    static constexpr size_t* hist  = nullptr;
};

template<size_t BS = 4, size_t NH = 3, size_t TL = 128,
         template <class> class DisStrat = dstrat_multi_bfs,
         class HistCount = hist_count>
struct Config
{
    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;
    static constexpr size_t nh = NH;

    template <class T>
    using DisStrat_temp = DisStrat<T>;

    using HistCount_t = HistCount;
};
