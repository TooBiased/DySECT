#pragma once

#include <cmath>
#include "cuckoo_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Config = CuckooConfig<> >
class TGrowingCuckoo : public CuckooTraits<TGrowingCuckoo<K,D,HF,Config> >::Base_t
{
private:
    using This_t         = TGrowingCuckoo<K,D,HF,Config>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = CuckooTraits<This_t>::Config_t::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::Config_t::tl;

    TGrowingCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                  size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        double avg_size_f = double(cap) * size_constraint / double(tl*bs);

        grow_amount    = 1;
        while(avg_size_f > (grow_amount << 1))
            grow_amount <<= 1;

        grow_table = (grow_amount < avg_size_f)
            ? std::floor(double(cap) * alpha / double(grow_amount * bs))-tl
            : 0;

        for (size_t i = grow_table; i < tl; ++i)
        {
            llb[i] = grow_amount - 1;
            llt[i] = std::make_unique<Bucket_t[]>(grow_amount);
        }

        for (size_t i = 0; i < grow_table; ++i)
        {
            llb[i] = (grow_amount << 1) - 1;
            llt[i] = std::make_unique<Bucket_t[]>(grow_amount << 1);
        }

        capacity    = (grow_table+tl) * grow_amount;
        grow_thresh = std::ceil((capacity + grow_amount*bs)/alpha);
    }

    TGrowingCuckoo(const TGrowingCuckoo&) = delete;
    TGrowingCuckoo& operator=(const TGrowingCuckoo&) = delete;

    TGrowingCuckoo(TGrowingCuckoo&& rhs)
        : Base_t(std::move(rhs)),
          grow_table(rhs.grow_table),
          grow_amount(rhs.grow_amount),
          grow_thresh(rhs.grow_thresh)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    TGrowingCuckoo& operator=(TGrowingCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        std::swap(grow_table , rhs.grow_table );
        std::swap(grow_amount, rhs.grow_amount);
        std::swap(grow_thresh, rhs.grow_thresh);

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(llb[i], rhs.llb[i]);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(llb[i]+1, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::alpha;
    using Base_t::capacity;

    using Base_t::h;

private:
    size_t                      llb[tl];
    std::unique_ptr<Bucket_t[]> llt[tl];

    size_t grow_table;
    size_t grow_amount;
    size_t grow_thresh;

    inline Bucket_t* getBucket1(HashSplitter_t h) const
    { return &(llt[h.tab1][h.loc1 & llb[h.tab1]]); }

    inline Bucket_t* getBucket2(HashSplitter_t h) const
    { return &(llt[h.tab2][h.loc2 & llb[h.tab2]]); }

    inline void inc_n()
    {
        ++n;
        if (n > grow_thresh) grow();
    }

    inline void grow()
    {
        size_t nsize = grow_amount << 1;
        auto   ntab  = std::make_unique<Bucket_t[]>(nsize);
        migrate(grow_table, ntab, nsize-1);

        llb[grow_table] = nsize-1;
        llt[grow_table] = std::move(ntab);

        capacity += grow_amount * bs;
        if (++grow_table == tl) { grow_table = 0; grow_amount <<= 1; }
        grow_thresh = std::ceil((capacity + grow_amount*bs)/alpha);
    }

    inline void migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, size_t bitmask)
    {
        for (size_t i = 0; i <= llb[tab]; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);
            for (size_t j = 0; j < bs; ++j)
            {
                auto e    = curr->elements[j];
                if (! e.first) break;
                auto hash = h(e.first);
                if      (getBucket1(hash) == curr) target[bitmask & hash.loc1].insert(e.first, e.second);
                else if (getBucket2(hash) == curr) target[bitmask & hash.loc2].insert(e.first, e.second);
                else
                {
                    std::cout << "something is wrong neither in first, nor second bucket." << std::endl;
                    exit(64);
                }
            }
        }
    }
};

template<class K, class D, class HF,
         class Config>
class CuckooTraits<TGrowingCuckoo<K,D,HF,Config> >
{
public:
    using Specialized_t  = TGrowingCuckoo<K,D,HF,Config>;
    using Base_t         = CuckooBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using HashFct_t      = HF;
    using Config_t       = Config;

    static constexpr size_t tl = Config::tl;
    static constexpr size_t bs = Config::bs;

    using Bucket_t       = Bucket<K,D,bs>;

    union HashSplitter_t
    {
        static constexpr size_t log(size_t k)
        { return (k-1) ? 1+log(k>>1) : 0; }

        static_assert( tl == 1ull<<log(tl),
                       "TL must be a power of two >0!");

        uint64_t hash;
        struct
        {
            uint64_t tab1 : log(tl);
            uint64_t tab2 : log(tl);
            uint64_t loc1 : 32-log(tl);
            uint64_t loc2 : 32-log(tl);
        };
    };
};


template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using GrowingCuckoo = TGrowingCuckoo<K,D,HF,CuckooConfig<BS,TL,DS,no_hist_count> >;

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using GrowingCuckooHist = TGrowingCuckoo<K,D,HF,CuckooConfig<BS,TL,DS,hist_count> >;
