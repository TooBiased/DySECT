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
    friend iterator_incr<This_t>;

private:
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;


    CuckooEG2L(size_t cap = 0      , double size_constraint = 1.1,
                       size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        double avg_size_f = double(cap) * size_constraint / double(tl*bs);

        size_t size_small    = 1;
        while(avg_size_f > (size_small << 1))
            size_small <<= 1;

        n_large = (size_small < avg_size_f)
            ? std::floor(double(cap) * alpha / double(size_small * bs))-tl
            : 0;

        for (size_t i = 0; i < n_large; ++i)
        {
            llt[i] = std::make_unique<Bucket_t[]>(size_small << 1);
        }

        for (size_t i = n_large; i < tl; ++i)
        {
            llt[i] = std::make_unique<Bucket_t[]>(size_small);
        }

        capacity    = (n_large+tl) * size_small * bs;
        bits_small  =  size_small - 1;
        bits_large  = (size_small << 1) - 1;
        grow_thresh = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh= 0; // ensures no shrinking until grown at least once
    }

    CuckooEG2L(const CuckooEG2L&) = delete;
    CuckooEG2L& operator=(const CuckooEG2L&) = delete;

    CuckooEG2L(CuckooEG2L&& rhs)
        : Base_t(std::move(rhs)),
          //grow_table(rhs.grow_table),
          //grow_amount(rhs.grow_amount),
          n_large(rhs.n_large),
          bits_small(rhs.bits_small),
          bits_large(rhs.bits_large),
          shrnk_thresh(rhs.shrnk_thresh)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            //llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    CuckooEG2L& operator=(CuckooEG2L&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        //std::swap(grow_table , rhs.grow_table );
        //std::swap(grow_amount, rhs.grow_amount);
        std::swap(n_large   , rhs.n_large);
        std::swap(bits_small, rhs.bits_small);
        std::swap(bits_large, rhs.bits_large);
        std::swap(shrnk_thresh,rhs.shrnk_thresh);

        for (size_t i = 0; i < tl; ++i)
        {
            //std::swap(llb[i], rhs.llb[i]);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(bitmask(i)+1, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;
    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;
    using Base_t::make_iterator;
    using Base_t::make_citerator;

    iterator begin()
    {
        auto temp = make_iterator(&llt[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    const_iterator begin() const
    {
        auto temp = make_citerator(&llt[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    std::unique_ptr<Bucket_t[]> llt[tl];

    size_t n_large;
    size_t bits_small;
    size_t bits_large;
    size_t shrnk_thresh;

    static constexpr size_t tl_bitmask = tl - 1;

    inline size_t bitmask(size_t tab) const
    { return (tab < n_large) ? bits_large : bits_small; }

    inline void      getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_t i) const
    {
        size_t tab = Ext::tab(h,i);
        size_t loc = Ext::loc(h,i) & bitmask(tab);
        return &(llt[tab][loc]);
    }

    inline void dec_n()
    {
        --n;
        if (n < shrnk_thresh) shrink();
    }

    inline void grow()
    {
        auto   ntab  = std::make_unique<Bucket_t[]>( bits_large + 1 );
        migrate_grw(n_large, ntab);

        llt[n_large] = std::move(ntab);

        capacity    += (bits_small+1) * bs;
        if (++n_large == tl)
        { n_large = 0; bits_small = bits_large; bits_large = (bits_large<<1) +1; }
        grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
    }

    inline void migrate_grw(size_t tab, std::unique_ptr<Bucket_t[]>& target)
    {
        size_t flag = bits_small+1;

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
                        (loc & bits_small) == i)
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
        if (n_large) { n_large--; }
        else         { n_large = tl-1; bits_small >>= 1; bits_large >>= 1; }
        auto ntab = std::make_unique<Bucket_t[]>(bits_small + 1);
        std::vector<std::pair<Key, Data> > buffer;

        migrate_shrnk( n_large, ntab, buffer );

        llt[n_large] = std::move(ntab);

        finish_shrnk(buffer);

        capacity    -= (bits_small+1) * bs;
        grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
        if (bits_small == 0 && !n_large) shrnk_thresh = 0;
    }

    inline void migrate_shrnk(size_t tab, std::unique_ptr<Bucket_t[]>& target,
                              std::vector<std::pair<Key, Data> >& buffer)
    {
        size_t flag = bits_small + 1;

        for (size_t i = 0; i < flag; ++i)
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
                        (Ext::loc(hash, ti) & bits_small) == i)
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
                             (Ext::loc(hash, ti) & bits_small) == i)
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
            bla += (Base_t::insert(e).second) ? 1: 0;
        }
    }

};



// Traits class defining types *************************************************

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



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooEG2L<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooEG2L<K,D,HF,Conf>;
    using ipointer = std::pair<const K,D>*;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t bs = Conf::bs;

public:
    iterator_incr(const Table_t& table_)
        : table(table_), end_tab(nullptr), tab(tl + 1)
    { }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    ipointer next(ipointer cur)
    {
        if (tab > tl) initialize_tab(cur);

        auto temp = cur+1;
        if (temp > end_tab)
        { temp = overflow_tab(); if (!temp) return nullptr; }

        while (!temp->first)
        {
            if (++temp > end_tab) return overflow_tab();
        }
        return temp;
    }

    size_t n_forwards;
private:
    const Table_t& table;
    ipointer       end_tab;
    size_t         tab;

    ipointer overflow_tab()
    {
        if (++tab >= tl) return nullptr;
        size_t size = (tab < table.n_large) ? table.bits_large : table.bits_small;
        end_tab = reinterpret_cast<ipointer>(&table.llt[tab][size].elements[bs-1]);
        return reinterpret_cast<ipointer>(&table.llt[tab][0].elements[0]);
    }

    void initialize_tab(ipointer ptr)
    {
        for (size_t i = 0; i < table.n_large; ++i)
        {
            ipointer tab_b_ptr = reinterpret_cast<ipointer>
                (&table.llt[i][0].elements[0]);
            ipointer tab_e_ptr = reinterpret_cast<ipointer>
                (&table.llt[i][table.bits_large].elements[bs-1]);

            if (tab_b_ptr <= ptr && ptr <= tab_e_ptr)
            {
                tab = i;
                end_tab = tab_e_ptr;
                return;
            }
        }
        for (size_t i = table.n_large; i < tl; ++i)
        {
            ipointer tab_b_ptr = reinterpret_cast<ipointer>
                (&table.llt[i][0].elements[0]);
            ipointer tab_e_ptr = reinterpret_cast<ipointer>
                (&table.llt[i][table.bits_small].elements[bs-1]);

            if (tab_b_ptr <= ptr && ptr <= tab_e_ptr)
            {
                tab = i;
                end_tab = tab_e_ptr;
                return;
            }
        }
    }

};
