#pragma once

/*******************************************************************************
 * include/cuckoo_eg2l.h
 *
 * Described in: https://arxiv.org/abs/1705.00997
 * Requirement:  OverAllocation (only for ceipg2l)
 *
 * Here we implement the main DySECT variants. Growability is
 * optained, by splitting the table into myltiple subtables thus
 * making it feasible to double the size of one subtable without
 * breaking a size constraint. Load imbalances between subtables are
 * efficiently reliefed by inter table exchanges.
 *
 * Two variants are implemented, eg2l is a variant that stores
 * independent pointers for each subtable migrations move elements to
 * a new table.  eipg2l grows each subtable in place, all subtables
 * are part of one large (overallocated chunk of memory).
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

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
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;
    using size_type      = typename Base_t::size_type;


    CuckooEG2L(size_type cap = 0      , double size_constraint = 1.1,
                       size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed)
    {
        double avg_size_f = double(cap) * size_constraint / double(tl*bs);

        size_type size_small    = 1;
        while(avg_size_f > (size_small << 1))
            size_small <<= 1;

        n_large = (size_small < avg_size_f)
            ? std::floor(double(cap) * alpha / double(size_small * bs))-tl
            : 0;

        for (size_type i = 0; i < n_large; ++i)
        {
            llt[i] = std::make_unique<Bucket_t[]>(size_small << 1);
        }

        for (size_type i = n_large; i < tl; ++i)
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
          n_large(rhs.n_large),
          bits_small(rhs.bits_small),
          bits_large(rhs.bits_large),
          shrnk_thresh(rhs.shrnk_thresh)
    {
        for (size_type i = 0; i < tl; ++i)
        {
            //llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    CuckooEG2L& operator=(CuckooEG2L&& rhs)
    {
        Base_t::operator=(std::move(rhs));
        std::swap(n_large   , rhs.n_large);
        std::swap(bits_small, rhs.bits_small);
        std::swap(bits_large, rhs.bits_large);
        std::swap(shrnk_thresh,rhs.shrnk_thresh);

        for (size_type i = 0; i < tl; ++i)
        {
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    static constexpr size_type bs = CuckooTraits<This_t>::bs;
    static constexpr size_type tl = CuckooTraits<This_t>::tl;
    static constexpr size_type nh = CuckooTraits<This_t>::nh;

    size_type n_large;
    size_type bits_small;
    size_type bits_large;
    size_type shrnk_thresh;

    std::unique_ptr<Bucket_t[]> llt[tl];

    static constexpr size_type tl_bitmask = tl - 1;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (i < tl) ? std::make_pair(bitmask(i)+1, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    iterator begin()
    {
        auto temp = make_iterator(&llt[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    const_iterator cbegin() const
    {
        auto temp = make_citerator(&llt[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    // Functions for finding buckets *******************************************

    inline size_type bitmask(size_type tab) const
    { return (tab < n_large) ? bits_large : bits_small; }

    inline void      getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_type i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_type i) const
    {
        size_type tab = Ext::tab(h,i);
        size_type loc = Ext::loc(h,i) & bitmask(tab);
        return &(llt[tab][loc]);
    }



    // Size changes (GROWING) **************************************************

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

    inline void migrate_grw(size_type tab, std::unique_ptr<Bucket_t[]>& target)
    {
        size_type flag = bits_small+1;

        for (size_type i = 0; i < flag; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);

            size_type tj0  = 0;
            Bucket_t* tar0 = &(target[i]);
            size_type tj1  = 0;
            Bucket_t* tar1 = &(target[i+flag]);

            for (size_type j = 0; j < bs; ++j)
            {
                auto e    = curr->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    size_type loc = Ext::loc(hash, ti);
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



    // Size changes (SHRINKING) ************************************************

    inline void dec_n()
    {
        --n;
        if (n < shrnk_thresh) shrink();
    }

    inline void shrink()
    {
        if (n_large) { n_large--; }
        else         { n_large = tl-1; bits_small >>= 1; bits_large >>= 1; }
        auto ntab = std::make_unique<Bucket_t[]>(bits_small + 1);
        std::vector<std::pair<key_type, mapped_type> > buffer;

        migrate_shrnk( n_large, ntab, buffer );

        llt[n_large] = std::move(ntab);

        finish_shrnk(buffer);

        capacity    -= (bits_small+1) * bs;
        grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
        if (bits_small == 0 && !n_large) shrnk_thresh = 0;
    }

    inline void migrate_shrnk(size_type tab, std::unique_ptr<Bucket_t[]>& target,
                              std::vector<std::pair<key_type, mapped_type> >& buffer)
    {
        size_type flag = bits_small + 1;

        for (size_type i = 0; i < flag; ++i)
        {
            Bucket_t* curr  = &(llt[tab][i]);
            Bucket_t* curr1 = &(llt[tab][i+flag]);
            Bucket_t* targ  = &(target[i]);
            size_type ind = 0;
            for (size_type j = 0; j < bs; ++j)
            {
                auto e = curr->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    if ( Ext::tab(hash, ti)  == tab &&
                        (Ext::loc(hash, ti) & bits_small) == i)
                    {
                        targ->elements[ind++] = e;
                        break;
                    }
                }
            }

            for (size_type j = 0; j < bs; ++j)
            {
                auto e = curr1->elements[j];
                if (! e.first)
                { break; }
                else if (ind >= bs)
                { buffer.push_back(e); }
                else
                {
                    auto hash = hasher(e.first);
                    for (size_type ti = 0; ti < nh; ++ti)
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

    inline void finish_shrnk(std::vector<std::pair<key_type, mapped_type> >& buffer)
    {
        size_type bla = 0;
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
    using Config_t       = Conf;

    using key_type       = K;
    using mapped_type    = D;
    using size_type      = size_t;

    static constexpr size_type tl = Config_t::tl;
    static constexpr size_type bs = Config_t::bs;
    static constexpr size_type nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooEG2L<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooEG2L<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using ipointer  = std::pair<const K,D>*;
    static constexpr size_type tl = Conf::tl;
    static constexpr size_type bs = Conf::bs;

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

private:
    const Table_t& table;
    ipointer       end_tab;
    size_type      tab;

    ipointer overflow_tab()
    {
        if (++tab >= tl) return nullptr;
        size_type size = (tab < table.n_large) ? table.bits_large : table.bits_small;
        end_tab = reinterpret_cast<ipointer>(&table.llt[tab][size].elements[bs-1]);
        return reinterpret_cast<ipointer>(&table.llt[tab][0].elements[0]);
    }

    void initialize_tab(ipointer ptr)
    {
        for (size_type i = 0; i < table.n_large; ++i)
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
        for (size_type i = table.n_large; i < tl; ++i)
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





// *****************************************************************************
// IN PLACE ********************************************************************
// *****************************************************************************

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooEIPG2L : public CuckooTraits<CuckooEIPG2L<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooEIPG2L<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;
    using size_type      = typename Base_t::size_type;

private:
    using value_intern   = std::pair<key_type, mapped_type>;

    static constexpr size_type bs = CuckooTraits<This_t>::bs;
    static constexpr size_type tl = CuckooTraits<This_t>::tl;
    static constexpr size_type nh = CuckooTraits<This_t>::nh;

    static constexpr size_type max_size     = 10ull << 30;
    static constexpr size_type max_loc_size = max_size / tl / sizeof(Bucket_t);
public:


    CuckooEIPG2L(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed)
    {
        auto temp = static_cast<Bucket_t*>(operator new (max_size));

        table     = std::unique_ptr<Bucket_t[]>(temp);

        double avg_size_f = double(cap) * size_constraint / double(tl*bs);

        size_type size_small = 1;
        while(avg_size_f > (size_small << 1))
            size_small <<= 1;

        n_large = (size_small < avg_size_f)
            ? std::floor(double(cap) * alpha / double(size_small * bs))-tl
            : 0;

        for (size_type i = 0; i < n_large; ++i)
        {
            // llt[i] = std::make_unique<Bucket_t[]>(size_small << 1);
            Bucket_t* offset = tableOff(i);
            std::fill(offset, offset + (size_small<<1), Bucket_t());
        }

        for (size_type i = n_large; i < tl; ++i)
        {
            // llt[i] = std::make_unique<Bucket_t[]>(size_small);
            auto offset = tableOff(i);
            std::fill(offset, offset + size_small, Bucket_t());
        }

        capacity    = (n_large+tl) * size_small * bs;
        bits_small  =  size_small - 1;
        bits_large  = (size_small << 1) - 1;
        grow_thresh = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh= 0; // ensures no shrinking until grown at least once
    }

    CuckooEIPG2L(const CuckooEIPG2L&)            = delete;
    CuckooEIPG2L& operator=(const CuckooEIPG2L&) = delete;

    CuckooEIPG2L(CuckooEIPG2L&&)                 = default;
    CuckooEIPG2L& operator=(CuckooEIPG2L&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    size_type n_large;
    size_type bits_small;
    size_type bits_large;
    size_type shrnk_thresh;

    std::unique_ptr<Bucket_t[]> table;

    static constexpr size_type tl_bitmask = tl - 1;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (i < tl) ? std::make_pair(bitmask(i)+1, tableOff(i))
                        : std::make_pair(0,nullptr);
    }

    iterator begin()
    {
        auto temp = make_iterator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    const_iterator cbegin() const
    {
        auto temp = make_citerator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    // Functions for finding buckets *******************************************

    inline size_type bitmask(size_type tab) const
    { return (tab < n_large) ? bits_large : bits_small; }

    inline void      getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_type i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_type i) const
    {
        size_type tab = Ext::tab(h,i);
        size_type loc = Ext::loc(h,i) & bitmask(tab);
        return tableOff(tab) + loc;
    }

    inline Bucket_t* tableOff(size_type t) const
    {
        return table.get() + t*max_loc_size;
    }

    // Size changes (GROWING) **************************************************

    void grow()
    {
        //auto   ntab  = std::make_unique<Bucket_t[]>( bits_large + 1 );
        Bucket_t* offset = tableOff(n_large);
        Bucket_t* new_s = offset + (bits_small + 1);
        Bucket_t* new_e = offset + (bits_large + 1);
        std::fill(new_s, new_e, Bucket_t());

        migrate_grw(n_large);
        //llt[n_large] = std::move(ntab);

        capacity    += (bits_small+1) * bs;
        if (++n_large == tl)
        { n_large = 0; bits_small = bits_large; bits_large = (bits_large<<1) +1; }
        grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
    }

    void migrate_grw(size_type tab)
    {
        size_type flag = bits_small+1;

        Bucket_t* b0 = tableOff(tab);
        Bucket_t* b1 = b0 + flag;

        for (size_type i = 0; i < flag; ++i, b0++, b1++)
        {
            size_type    k0  = 0;
            size_type    k1  = 0;

            for (size_type j = 0; j < bs; ++j)
            {
                auto e    = b0->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    size_type loc = Ext::loc(hash, ti);
                    if ( Ext::tab(hash, ti) == tab &&
                        (loc & bits_small) == i)
                    {
                        if (loc & flag) b1->elements[k1++] = e;
                        else            b0->elements[k0++] = e;
                        break;
                    }
                }
            }
            for (size_type j = k0; j < bs; ++j)
            {
                b0->elements[j] = value_intern();
            }
        }
    }



    // Size changes (SHRINKING) ************************************************

    // inline void dec_n()
    // {
    //     --n;
    //     if (n < shrnk_thresh) shrink();
    // }

    // inline void shrink()
    // {
    //     if (n_large) { n_large--; }
    //     else         { n_large = tl-1; bits_small >>= 1; bits_large >>= 1; }
    //     auto ntab = std::make_unique<Bucket_t[]>(bits_small + 1);
    //     std::vector<std::pair<key_type, mapped_type> > buffer;

    //     migrate_shrnk( n_large, ntab, buffer );

    //     llt[n_large] = std::move(ntab);

    //     finish_shrnk(buffer);

    //     capacity    -= (bits_small+1) * bs;
    //     grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
    //     shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
    //     if (bits_small == 0 && !n_large) shrnk_thresh = 0;
    // }

    // inline void migrate_shrnk(size_type tab, std::unique_ptr<Bucket_t[]>& target,
    //                           std::vector<std::pair<key_type, mapped_type> >& buffer)
    // {
    //     size_type flag = bits_small + 1;

    //     for (size_type i = 0; i < flag; ++i)
    //     {
    //         Bucket_t* curr  = &(llt[tab][i]);
    //         Bucket_t* curr1 = &(llt[tab][i+flag]);
    //         Bucket_t* targ  = &(target[i]);
    //         size_type ind = 0;
    //         for (size_type j = 0; j < bs; ++j)
    //         {
    //             auto e = curr->elements[j];
    //             if (! e.first) break;
    //             auto hash = hasher(e.first);

    //             for (size_type ti = 0; ti < nh; ++ti)
    //             {
    //                 if ( Ext::tab(hash, ti)  == tab &&
    //                     (Ext::loc(hash, ti) & bits_small) == i)
    //                 {
    //                     targ->elements[ind++] = e;
    //                     break;
    //                 }
    //             }
    //         }

    //         for (size_type j = 0; j < bs; ++j)
    //         {
    //             auto e = curr1->elements[j];
    //             if (! e.first)
    //             { break; }
    //             else if (ind >= bs)
    //             { buffer.push_back(e); }
    //             else
    //             {
    //                 auto hash = hasher(e.first);
    //                 for (size_type ti = 0; ti < nh; ++ti)
    //                 {
    //                     if ( Ext::tab(hash, ti)  == tab &&
    //                          (Ext::loc(hash, ti) & bits_small) == i)
    //                     {
    //                         targ->elements[ind++] = e;
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // inline void finish_shrnk(std::vector<std::pair<key_type, mapped_type> >& buffer)
    // {
    //     size_type bla = 0;
    //     n -= buffer.size();
    //     for (auto& e : buffer)
    //     {
    //         bla += (Base_t::insert(e).second) ? 1: 0;
    //     }
    // }

};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooEIPG2L<K,D,HF,Conf> >
{
public:
    using Specialized_t  = CuckooEIPG2L<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Config_t       = Conf;

    using key_type       = K;
    using mapped_type    = D;
    using size_type      = size_t;

    static constexpr size_type tl = Config_t::tl;
    static constexpr size_type bs = Config_t::bs;
    static constexpr size_type nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooEIPG2L<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooEIPG2L<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using ipointer  = std::pair<const K,D>*;
    static constexpr size_type tl = Conf::tl;
    static constexpr size_type bs = Conf::bs;

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

private:
    const Table_t& table;
    ipointer       end_tab;
    size_type      tab;

    ipointer overflow_tab()
    {
        if (++tab >= tl) return nullptr;
        auto tab_start = table.table.get() + tab*table.max_loc_size;
        size_type size = (tab < table.n_large) ? table.bits_large : table.bits_small;
        end_tab = reinterpret_cast<ipointer>(&(tab_start + size)->elements[bs-1]);
        return reinterpret_cast<ipointer>(&tab_start->elements[0]);
    }

    void initialize_tab(ipointer ptr)
    {
        auto tab_curr = table.table.get();
        if (ptr < reinterpret_cast<ipointer>(&tab_curr->elements[0])) return;

        for (size_type i = 0; i < tl; ++i)
        {
            auto tab_size = (i < table.n_large) ? table.bits_large : table.bits_small;
            auto bkt_last = tab_curr + tab_size;
            auto par_last = reinterpret_cast<ipointer>(& bkt_last->elements[bs-1]);

            //tab_offset += table.max_loc_size;
            if (ptr <= par_last)
            {
                tab = i;
                end_tab = par_last;
                return;
            }
        }
    }

};
