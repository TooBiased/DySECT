#pragma once

#include <cmath>
#include "cuckoo_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooEG2L : public CuckooTraits<CuckooEG2L<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooEG2L<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    friend Base_t;

public:
    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;

private:
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    CuckooEG2L(size_t cap = 0      , double size_constraint = 1.1,
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
        grow_thresh = std::ceil((capacity + 2*grow_amount*bs)/alpha);
        shrnk_thresh= 0; // ensures no shrinking until grown at least once
    }

    CuckooEG2L(const CuckooEG2L&) = delete;
    CuckooEG2L& operator=(const CuckooEG2L&) = delete;

    CuckooEG2L(CuckooEG2L&& rhs)
        : Base_t(std::move(rhs)),
          grow_table(rhs.grow_table),
          grow_amount(rhs.grow_amount),
          grow_thresh(rhs.grow_thresh),
          shrnk_thresh(rhs.shrnk_thresh)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    CuckooEG2L& operator=(CuckooEG2L&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        std::swap(grow_table , rhs.grow_table );
        std::swap(grow_amount, rhs.grow_amount);
        std::swap(grow_thresh, rhs.grow_thresh);
        std::swap(shrnk_thresh,rhs.shrnk_thresh);

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

    using Base_t::hasher;

private:
    size_t                      llb[tl];
    std::unique_ptr<Bucket_t[]> llt[tl];

    size_t grow_table;
    size_t grow_amount;
    size_t grow_thresh;
    size_t shrnk_thresh;


    static constexpr size_t tl_bitmask = tl - 1;

    inline void      getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_t i) const
    {
        size_t tab = Ext::tab(h,i);
        size_t loc = Ext::loc(h,i) & llb[tab];
        return &(llt[tab][loc]);
    }

    inline void inc_n()
    {
        ++n;
        if (n > grow_thresh) grow();
    }

    inline void dec_n()
    {
        --n;
        if (n < shrnk_thresh) shrink();
    }

    inline void grow()
    {
        size_t nsize = grow_amount << 1;
        auto   ntab  = std::make_unique<Bucket_t[]>(nsize);
        migrate_grw(grow_table, ntab);

        llb[grow_table] = nsize-1;
        llt[grow_table] = std::move(ntab);

        capacity    += grow_amount * bs;
        if (++grow_table == tl) { grow_table = 0; grow_amount <<= 1; }
        grow_thresh  = std::ceil((capacity + 2*grow_amount*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - 2*grow_amount*bs)/alpha);
    }

    inline void migrate_grw(size_t tab, std::unique_ptr<Bucket_t[]>& target)
    {
        size_t bold = llb[tab];
        size_t flag = bold+1;

        for (size_t i = 0; i < flag; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);

            size_t    tj0  = 0;
            Bucket_t* tar0 = &(target[i]);
            size_t    tj1  = 0;
            Bucket_t* tar1 = &(target[i+flag]);

            for (size_t j = 0; j < bs; ++j)
            {
                auto e    = curr->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    size_t loc = Ext::loc(hash, ti);
                    if ( Ext::tab(hash, ti) == tab &&
                        (loc & bold) == i)
                    {
                        if (loc & flag) tar1->elements[tj1++] = e;
                        else            tar0->elements[tj0++] = e;
                        break;
                    }
                }
            }
        }
    }

    inline void shrink()
    {
        if (grow_table) { grow_table--; }
        else            { grow_table = tl-1; grow_amount >>= 1; }
        auto ntab = std::make_unique<Bucket_t[]>(grow_amount);
        std::vector<std::pair<Key, Data> > buffer;

        migrate_shrnk( grow_table, ntab, buffer );

        llb[grow_table] = grow_amount-1;
        llt[grow_table] = std::move(ntab);

        finish_shrnk(buffer);

        capacity    -= grow_amount * bs;
        grow_thresh  = std::ceil((capacity + 2*grow_amount*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - 2*grow_amount*bs)/alpha);
        if (grow_amount == 1 && !grow_table) shrnk_thresh = 0;
    }

    inline void migrate_shrnk(size_t tab, std::unique_ptr<Bucket_t[]>& target,
                              std::vector<std::pair<Key, Data> >& buffer)
    {
        size_t bnew = llb[tab] >> 1;
        size_t flag = bnew + 1;

        for (size_t i = 0; i <= bnew; ++i)
        {
            Bucket_t* curr  = &(llt[tab][i]);
            Bucket_t* curr1 = &(llt[tab][i+flag]);
            Bucket_t* targ  = &(target[i]);
            size_t ind = 0;
            for (size_t j = 0; j < bs; ++j)
            {
                auto e = curr->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if ( Ext::tab(hash, ti)  == tab &&
                        (Ext::loc(hash, ti) & bnew) == i)
                    {
                        targ->elements[ind++] = e;
                        break;
                    }
                }
            }

            for (size_t j = 0; j < bs; ++j)
            {
                auto e = curr1->elements[j];
                if (! e.first)
                { break; }
                else if (ind >= bs)
                { buffer.push_back(e); }
                else
                {
                    auto hash = hasher(e.first);
                    for (size_t ti = 0; ti < nh; ++ti)
                    {
                        if ( Ext::tab(hash, ti)  == tab &&
                             (Ext::loc(hash, ti) & bnew) == i)
                        {
                            targ->elements[ind++] = e;
                            break;
                        }
                    }
                }
            }

        }
    }

    inline void finish_shrnk(std::vector<std::pair<Key, Data> >& buffer)
    {
        size_t bla = 0;
        n -= buffer.size();
        for (auto& e : buffer)
        {
            bla += (Base_t::insert(e)) ? 1: 0;
        }
    }

};

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooEG2L<K,D,HF,Conf> >
{
public:
    using Specialized_t  = CuckooEG2L<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Config_t       = Conf;

    static constexpr size_t tl = Config_t::tl;
    static constexpr size_t bs = Config_t::bs;
    static constexpr size_t nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;
};
