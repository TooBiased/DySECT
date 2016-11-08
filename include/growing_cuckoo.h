#pragma once

#include <cmath>
#include "cuckoo_base.h"
#include "strategies/dstrat_triv.h"

template<size_t FL>
union TwoLvlHashSplitter
{
    static constexpr size_t log(size_t k)
    {
        return (k-1) ? 1+log(k>>1) : 0;
    }
    static_assert( FL == 1ull<<log(FL),
                   "FL must be a power of two >0!");

    uint64_t hash;
    struct
    {
        uint64_t tab1 : log(FL);
        uint64_t tab2 : log(FL);
        uint64_t loc1 : 32-log(FL);
        uint64_t loc2 : 32-log(FL);
    };
};

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
class GrowingCuckoo : public CuckooBase<K, D, HF, DS, BS,
                                        GrowingCuckoo<K,D,HF,DS,TL,BS>,
                                        TwoLvlHashSplitter<TL> >
{
private:
    using This_t       = GrowingCuckoo<K,D,HF,DS,TL,BS>;
    using Base_t       = CuckooBase<K,D,HF,DS,BS,This_t,TwoLvlHashSplitter<TL> >;
    friend Base_t;
    using Bucket_t     = typename Base_t::Bucket_t;
    using HashSplitter = typename Base_t::HashSplitter;

public:
    GrowingCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                  size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        std::cout << cap << std::endl;

        double avg_size_f = double(cap) * size_constraint / double(TL*BS);
        std::cout << avg_size_f << std::endl;

        grow_amount    = 1;
        while(avg_size_f > (grow_amount << 1))
            grow_amount <<= 1;

        std::cout << grow_amount << std::endl;

        grow_table = (grow_amount < avg_size_f) ? std::floor(double(cap) * alpha / double(grow_amount * BS))-TL : 0;

        std::cout << grow_table << std::endl;

        for (size_t i = grow_table; i < TL; ++i)
        {
            llb[i] = grow_amount - 1;
            llt[i] = std::make_unique<Bucket_t[]>(grow_amount);
        }

        for (size_t i = 0; i < grow_table; ++i)
        {
            llb[i] = (grow_amount << 1) - 1;
            llt[i] = std::make_unique<Bucket_t[]>(grow_amount << 1);
        }

        capacity    = (grow_table+TL) * grow_amount;
        grow_thresh = std::ceil((capacity + grow_amount*BS)/alpha);
    }

    GrowingCuckoo(const GrowingCuckoo&) = delete;
    GrowingCuckoo& operator=(const GrowingCuckoo&) = delete;

    GrowingCuckoo(GrowingCuckoo&& rhs)
        : Base_t(std::move(rhs)),
          grow_table(rhs.grow_table),
          grow_amount(rhs.grow_amount),
          grow_thresh(rhs.grow_thresh)
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    GrowingCuckoo& operator=(GrowingCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        std::swap(grow_table , rhs.grow_table );
        std::swap(grow_amount, rhs.grow_amount);
        std::swap(grow_thresh, rhs.grow_thresh);

        for (size_t i = 0; i < TL; ++i)
        {
            std::swap(llb[i], rhs.llb[i]);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < TL) ? std::make_pair(llb[i]+1, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    static constexpr size_t tl = TL;

    using Base_t::n;
    using Base_t::alpha;
    using Base_t::capacity;

    using Base_t::h;

private:
    size_t                      llb[TL];
    std::unique_ptr<Bucket_t[]> llt[TL];

    size_t grow_table;
    size_t grow_amount;
    size_t grow_thresh;

    inline Bucket_t* getBucket1(HashSplitter h) const
    { return &(llt[h.tab1][h.loc1 & llb[h.tab1]]); }

    inline Bucket_t* getBucket2(HashSplitter h) const
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

        capacity += grow_amount * BS;
        if (++grow_table == TL) { grow_table = 0; grow_amount <<= 1; }
        grow_thresh = std::ceil((capacity + grow_amount*BS)/alpha);
    }

    inline void migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, size_t bitmask)
    {
        for (size_t i = 0; i <= llb[tab]; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);
            for (size_t j = 0; j < BS; ++j)
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
