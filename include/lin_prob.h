#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "cuckoo_base.h"
/*
class hist_count
{
public:
    hist_count(size_t s) : steps(s), hist(new size_t[s])
    { for (size_t i = 0; i < s; ++i) { hist[i] = 0; } }

    void add(size_t i) { ++hist[i];}

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


template<size_t BS = 8, size_t TL = 256,
         template <class> class DisStrat = dstrat_triv,
         class HistCount = no_hist_count>
struct CuckooConfig
{
    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;

    template <class T>
    using DisStrat_temp = DisStrat<T>;

    using HistCount_t = HistCount;
};

template<class T>
class CuckooTraits;
 EXAMPLE IMPLEMENTATION
{
public:
    using Specialized_t  = T;
    using Base_t         = CuckooBase<T>;
    using Key            = ... ;
    using Data           = ... ;
    using Bucket_t       = Bucket<Key,Data,bs>;
    using HashFct_t      = ... ;

    using Config_t       = CuckooConfig<...>;

    union HashSplitter
    {
        uint64_t hash;
        struct
        {
            ... //partial hash bits;
        };
    };
};*/


template<class K, class D, class HF = std::hash<K>, class Config = CuckooConfig<> >
class FastLinProbTable
{
private:

    using This_t         = FastLinProbTable<K,D,HF,Config>;
    using HashFct_t      = HF;
    using HistCount_t    = typename Config::HistCount_t;
    using Bucket_t       = Bucket<K,D,1>;
    using Cell_t         = std::pair<K,D>;

public:
    using Key      = K;
    using Data     = D;
    using FRet     = std::pair<bool, Data>;

    FastLinProbTable(size_t cap = 0, double size_constraint = 1.1,
                     size_t dis_steps = 0, size_t /*unneeded*/ = 0)
        : n(0), alpha(size_constraint), steps(dis_steps), hcounter(dis_steps)
    {
        capacity = Config::tl;
        while( capacity < cap ) capacity <<= 1;
        capacity <<= 1;
        table = std::make_unique<Cell_t[]>(capacity);
        grow_thresh = capacity * 0.6;
        bitmask = capacity - 1;
    }

    ~FastLinProbTable() = default;

    FastLinProbTable(const FastLinProbTable&) = delete;
    FastLinProbTable& operator=(const FastLinProbTable&) = delete;

    FastLinProbTable(FastLinProbTable&& rhs)
        : n(rhs.n), capacity(rhs.capacity), bitmask(rhs.bitmask), alpha(rhs.alpha),
          steps(rhs.steps), grow_thresh(rhs.grow_thresh),
          hcounter(rhs.steps), table(std::move(rhs.table))
    { }

    FastLinProbTable& operator=(FastLinProbTable&& rhs)
    {
        n = rhs.n;   capacity = rhs.capacity;   bitmask = rhs.bitmask;
        alpha = rhs.alpha;  grow_thresh = rhs.grow_thresh;
        table = std::move(rhs.table);
        return *this;
    }

    bool insert(Key k, Data d);
    bool insert(std::pair<Key,Data> t);
    FRet find  (Key k) const;
    bool remove(Key k);

    /*** some functions for easier load visualization (not in final product) **/
    std::pair<size_t, Bucket_t*> getTable(size_t)
    { return std::make_pair(0, nullptr); }
    void clearHist()
    { for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0; }

    /*** members that should become private at some point *********************/
    size_t     n;
    size_t     capacity;
    size_t     bitmask;
    double     alpha;
    size_t     steps;
    size_t     grow_thresh;
    HashFct_t  hasher;
    HistCount_t  hcounter;
    std::unique_ptr<Cell_t[]> table;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

private:
    /*** static polymorph functions *******************************************/
    inline void inc_n() { ++n; if (n > grow_thresh) grow(); }
    inline void dec_n() { --n; }

    inline size_t hash(Key k) const
    { return hasher(k) & bitmask; }

    inline bool before(size_t i0, size_t i1)
    {
        if (i1 < i0) return false;
        if (i0 < i1 && i1 < bitmask) return true;
    }

    FastLinProbTable(size_t cap, size_t steps)
        : n(0), capacity(cap), bitmask(cap-1), alpha(0), steps(steps),
          grow_thresh(cap), hcounter(steps), table(new Cell_t[cap])
    { }

    inline void grow()
    {
        auto nsize  = (bitmask+1) << 1;
        auto ntable = This_t(nsize, steps);

        for (size_t i = 0; i <= bitmask; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        std::swap(bitmask, ntable.bitmask);
        std::swap(table  , ntable.table);
        grow_thresh = (bitmask+1) * 0.6;
    }

};



/* IMPLEMENTATION *************************************************************/
template<class K, class D, class HF, class Config>
inline bool FastLinProbTable<K,D,HF,Config>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class K, class D, class HF, class Config>
inline bool FastLinProbTable<K,D,HF,Config>::insert(std::pair<Key, Data> t)
{
    auto ind = hash(t.first);

    for (size_t i = ind; i < ind+steps; ++i)
    {
        auto temp = table[i&bitmask];
        if ( temp.first == 0 )
        {
            table[i&bitmask] = t;
            hcounter.add(i - ind);
            inc_n();
            return true;
        }
    }
    return false;
}

template<class K, class D, class HF, class Config>
inline typename FastLinProbTable<K,D,HF,Config>::FRet
FastLinProbTable<K,D,HF,Config>::find(Key k) const
{
    auto ind = hash(k);

    for (size_t i = ind; i < ind+steps; ++i)
    {
        auto temp = table[i & bitmask];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            return std::make_pair(true, temp.second);
        }
    }
    return std::make_pair(false, 0);
}

template<class K, class D, class HF, class Config>
inline bool FastLinProbTable<K,D,HF,Config>::remove(Key k)
{
    auto ind = hash(k);
    size_t i = ind;
    while (true)
    {
        auto temp = table[i];
        if ( temp.first == 0 )
        {
            return false;
        }
        else if ( temp.first == k )
        {
            auto i2 = i+1;
            while (true)
            {
                auto temp2 = table[i2 & bitmask];
                if (temp2.first == 0)
                    break;
                else if (before(hash(temp2.first), i))
                {
                    table[i & bitmask] = temp2;
                    i                  = i2;
                }
                i2++;
            }
            table[i & bitmask] = Cell_t(0, 0);
        }

        ++i;
    }
}
