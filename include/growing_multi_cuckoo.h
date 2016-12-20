#pragma once

#include <cmath>
#include "cuckoo_multi_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class GrowingMultiCuckoo : public CuckooTraits<GrowingMultiCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = GrowingMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    friend Base_t;

    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;

    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    //using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = Hasher_t::Hashed_t;
    using Extractor_t    = Hasher_t::template Extractor<nh>;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    GrowingMultiCuckoo(size_t cap = 0      , double size_constraint = 1.1,
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

        capacity    = (grow_table+tl) * grow_amount * bs;
        grow_thresh = std::ceil((capacity + grow_amount*bs)/alpha);
    }

    GrowingMultiCuckoo(const GrowingMultiCuckoo&) = delete;
    GrowingMultiCuckoo& operator=(const GrowingMultiCuckoo&) = delete;

    GrowingMultiCuckoo(GrowingMultiCuckoo&& rhs)
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

    GrowingMultiCuckoo& operator=(GrowingMultiCuckoo&& rhs)
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


    static constexpr size_t tl_bitmask = tl - 1;

    inline void      getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_t i) const
    {
        size_t tab = Extractor_t::tab(h,i); //(h.tab1+i*h.tab2) & tl_bitmask;
        size_t loc = Extractor_t::loc(h,i) & llb[tab]; //(h.loc1+i*h.loc2) & llb[tab];
        return &(llt[tab][loc]);
    }

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
                auto hash = hasher(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (  Extractor_t::tab(hash, ti) == tab && //hash.tab1+ti*hash.tab2 == tab &&
                         (Extractor_t::loc(hash, ti) & llb[tab])  == i) //((hash.loc1+ti*hash.loc2) & llb[tab]) == i)
                    {
                        target[Extractor_t::loc(hash, ti) & bitmask].insert(e.first, e.second);
                        break;
                    }
                }
            }
        }
    }
};

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<GrowingMultiCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = GrowingMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Config_t       = Conf;
    //using HashFct_t      = HF;

    static constexpr size_t tl = Config_t::tl;
    static constexpr size_t bs = Config_t::bs;
    static constexpr size_t nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), 32-ct_log(tl), 2, 1>;
    using Bucket_t       = Bucket<K,D,bs>;

    /*union HashSplitter_t
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
    */
};
