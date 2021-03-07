#pragma once

/*******************************************************************************
 * include/cuckoo_dysect.h
 *
 * Described in: https://arxiv.org/abs/1705.00997
 * Requirement:  OverAllocation (only for dysect_inplace)
 *
 * Here we implement the main DySECT variants. Growability is
 * optained, by splitting the table into myltiple subtables thus
 * making it feasible to double the size of one subtable without
 * breaking a size constraint. Load imbalances between subtables are
 * efficiently reliefed by inter table exchanges.
 *
 * Two variants are implemented, cuckoo_dysect is a variant that
 * stores independent pointers for each subtable migrations move
 * elements to a new table.  cuckoo_dysect_inplace grows each subtable
 * in place, all subtables are part of one large (overallocated chunk
 * of memory).  This can be more efficient since offsets can be
 * computed more quickly.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <cmath>
#include "utils/default_hash.hpp"
#include "cuckoo_base.h"

namespace dysect
{
    template<class T>
    class cuckoo_traits;

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_dysect : public cuckoo_traits<cuckoo_dysect<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type      = cuckoo_dysect<K,D,HF,Conf>;
        using base_type      = typename cuckoo_traits<this_type>::base_type;
        using bucket_type    = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type    = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type    = typename hasher_type::hashed_type;
        using ext            = typename hasher_type::extractor_type;

        friend base_type;
        friend iterator_incr<this_type>;

    public:
        using key_type       = typename cuckoo_traits<this_type>::key_type;
        using mapped_type    = typename cuckoo_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;
        using size_type      = typename base_type::size_type;


        cuckoo_dysect(size_type cap = 0      , double size_constraint = 1.1,
                      size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed)
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
                llt[i] = std::make_unique<bucket_type[]>(size_small << 1);
            }

            for (size_type i = n_large; i < tl; ++i)
            {
                llt[i] = std::make_unique<bucket_type[]>(size_small);
            }

            capacity    = (n_large+tl) * size_small * bs;
            bits_small  =  size_small - 1;
            bits_large  = (size_small << 1) - 1;

            if (n_large == tl)
            { n_large = 0; bits_small = bits_large; bits_large = (bits_large<<1) +1; }

            grow_thresh = std::ceil((capacity + (bits_large+1)*bs)/alpha);
            shrnk_thresh= 0; // ensures no shrinking until grown at least once
        }

        cuckoo_dysect(const cuckoo_dysect&) = delete;
        cuckoo_dysect& operator=(const cuckoo_dysect&) = delete;

        cuckoo_dysect(cuckoo_dysect&& rhs)
            : base_type(std::move(rhs)),
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

        cuckoo_dysect& operator=(cuckoo_dysect&& rhs)
        {
            base_type::operator=(std::move(rhs));
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
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        static constexpr size_type bs = cuckoo_traits<this_type>::bs;
        static constexpr size_type tl = cuckoo_traits<this_type>::tl;
        static constexpr size_type nh = cuckoo_traits<this_type>::nh;

        size_type n_large;
        size_type bits_small;
        size_type bits_large;
        size_type shrnk_thresh;

        std::unique_ptr<bucket_type[]> llt[tl];

        static constexpr size_type tl_bitmask = tl - 1;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
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

        inline size_type    bitmask(size_type tab) const
        { return (tab < n_large) ? bits_large : bits_small; }

        inline void         get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
                mem[i] = get_bucket(h,i);
        }

        inline bucket_type* get_bucket (hashed_type h, size_type i) const
        {
            size_type tab = ext::tab(h,i);
            size_type loc = ext::loc(h,i) & bitmask(tab);
            return &(llt[tab][loc]);
        }



        // Size changes (GROWING) **************************************************

        inline void grow()
        {
            auto   ntab  = std::make_unique<bucket_type[]>( bits_large + 1 );
            migrate_grw(n_large, ntab);

            llt[n_large] = std::move(ntab);

            capacity    += (bits_small+1) * bs;
            if (++n_large == tl)
            { n_large = 0; bits_small = bits_large; bits_large = (bits_large<<1) +1; }
            grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
            shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
        }

        inline void migrate_grw(size_type tab, std::unique_ptr<bucket_type[]>& target)
        {
            size_type flag = bits_small+1;

            for (size_type i = 0; i < flag; ++i)
            {
                bucket_type* curr = &(llt[tab][i]);

                size_type tj0  = 0;
                bucket_type* tar0 = &(target[i]);
                size_type tj1  = 0;
                bucket_type* tar1 = &(target[i+flag]);

                for (size_type j = 0; j < bs; ++j)
                {
                    auto e    = curr->elements[j];
                    if (! e.first) break;
                    auto hash = hasher(e.first);

                    for (size_type ti = 0; ti < nh; ++ti)
                    {
                        size_type loc = ext::loc(hash, ti);
                        if ( ext::tab(hash, ti) == tab &&
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
            auto ntab = std::make_unique<bucket_type[]>(bits_small + 1);
            std::vector<std::pair<key_type, mapped_type> > buffer;

            migrate_shrnk( n_large, ntab, buffer );

            llt[n_large] = std::move(ntab);

            finish_shrnk(buffer);

            capacity    -= (bits_small+1) * bs;
            grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
            shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
            if (bits_small == 0 && !n_large) shrnk_thresh = 0;
        }

        inline void migrate_shrnk(size_type tab, std::unique_ptr<bucket_type[]>& target,
                                  std::vector<std::pair<key_type, mapped_type> >& buffer)
        {
            size_type flag = bits_small + 1;

            for (size_type i = 0; i < flag; ++i)
            {
                bucket_type* curr  = &(llt[tab][i]);
                bucket_type* curr1 = &(llt[tab][i+flag]);
                bucket_type* targ  = &(target[i]);
                size_type ind = 0;
                for (size_type j = 0; j < bs; ++j)
                {
                    auto e = curr->elements[j];
                    if (! e.first) break;
                    auto hash = hasher(e.first);

                    for (size_type ti = 0; ti < nh; ++ti)
                    {
                        if ( ext::tab(hash, ti)  == tab &&
                             (ext::loc(hash, ti) & bits_small) == i)
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
                            if ( ext::tab(hash, ti)  == tab &&
                                 (ext::loc(hash, ti) & bits_small) == i)
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
            size_type err = 0;
            n -= buffer.size();
            for (auto& e : buffer)
            {
                err += (base_type::insert(e).second) ? 1: 0;
            }
        }

    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_dysect<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_dysect<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using key_type         = K;
        using mapped_type      = D;
        using size_type        = size_t;

        static constexpr size_type tl = config_type::tl;
        static constexpr size_type bs = config_type::bs;
        static constexpr size_type nh = config_type::nh;
        static constexpr bool fix_errors = config_type::fix_errors;

        using hasher_type    = hasher<K, HF, ct_log(tl), nh, true, true>;
        using bucket_type    = bucket<K,D,bs>;

    };


// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_dysect<K,D,HF,Conf> >
    {
    public:
        using table_type = cuckoo_dysect<K,D,HF,Conf>;
    private:
        using size_type  = typename table_type::size_type;
        using ipointer   = std::pair<const K,D>*;
        static constexpr size_type tl = Conf::tl;
        static constexpr size_type bs = Conf::bs;

    public:

        iterator_incr(const table_type& table_)
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
        const table_type& table;
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

    #include <stdlib.h>

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_dysect_inplace : public cuckoo_traits<cuckoo_dysect_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type      = cuckoo_dysect_inplace<K,D,HF,Conf>;
        using base_type      = typename cuckoo_traits<this_type>::base_type;
        using bucket_type    = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type    = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type    = typename hasher_type::hashed_type;
        using ext            = typename hasher_type::extractor_type;

        friend base_type;
        friend iterator_incr<this_type>;

    public:
        using key_type       = typename cuckoo_traits<this_type>::key_type;
        using mapped_type    = typename cuckoo_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;
        using size_type      = typename base_type::size_type;

    private:
        using value_intern   = std::pair<key_type, mapped_type>;

        static constexpr size_type bs = cuckoo_traits<this_type>::bs;
        static constexpr size_type tl = cuckoo_traits<this_type>::tl;
        static constexpr size_type nh = cuckoo_traits<this_type>::nh;

        static constexpr size_type max_size     = 10ull << 30;
        static constexpr size_type max_loc_size = max_size / tl / sizeof(bucket_type);
    public:


        cuckoo_dysect_inplace(size_type cap = 0      , double size_constraint = 1.1,
                              size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed)
        {
            //auto temp = operator new (max_size);
            auto temp = static_cast<bucket_type*>(aligned_alloc(4096,max_size));
            table     = std::unique_ptr<bucket_type[]>(temp);

            double avg_size_f = double(cap) * size_constraint / double(tl*bs);

            size_type size_small = 1;
            while(avg_size_f > (size_small << 1))
                size_small <<= 1;

            n_large = (size_small < avg_size_f)
                ? std::floor(double(cap) * alpha / double(size_small * bs))-tl
                : 0;

            for (size_type i = 0; i < n_large; ++i)
            {
                // llt[i] = std::make_unique<bucket_type[]>(size_small << 1);
                bucket_type* offset = table_off(i);
                std::fill(offset, offset + (size_small<<1), bucket_type());
            }

            for (size_type i = n_large; i < tl; ++i)
            {
                // llt[i] = std::make_unique<bucket_type[]>(size_small);
                auto offset = table_off(i);
                std::fill(offset, offset + size_small, bucket_type());
            }

            capacity    = (n_large+tl) * size_small * bs;
            bits_small  =  size_small - 1;
            bits_large  = (size_small << 1) - 1;

            if (n_large == tl)
            { n_large = 0; bits_small = bits_large; bits_large = (bits_large<<1) +1; }

            grow_thresh = std::ceil((capacity + (bits_large+1)*bs)/alpha);
            shrnk_thresh= 0; // ensures no shrinking until grown at least once
        }

        cuckoo_dysect_inplace(const cuckoo_dysect_inplace&)            = delete;
        cuckoo_dysect_inplace& operator=(const cuckoo_dysect_inplace&) = delete;

        cuckoo_dysect_inplace(cuckoo_dysect_inplace&&)                 = default;
        cuckoo_dysect_inplace& operator=(cuckoo_dysect_inplace&&) = default;

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        size_type n_large;
        size_type bits_small;
        size_type bits_large;
        size_type shrnk_thresh;

        //std::unique_ptr<bucket_type[]> memory;
        //bucket_type* table;
        std::unique_ptr<bucket_type[]> table;

        static constexpr size_type tl_bitmask = tl - 1;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
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

        inline size_type    bitmask(size_type tab) const
        { return (tab < n_large) ? bits_large : bits_small; }

        inline void         get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
                mem[i] = get_bucket(h,i);
        }

        inline bucket_type* get_bucket (hashed_type h, size_type i) const
        {
            size_type tab = ext::tab(h,i);
            size_type loc = ext::loc(h,i) & bitmask(tab);
            return table_off(tab) + loc;
        }

        inline bucket_type* table_off(size_type t) const
        {
            return table.get() + t*max_loc_size;
        }

        // Size changes (GROWING) **************************************************

        void grow()
        {
            //auto   ntab  = std::make_unique<bucket_type[]>( bits_large + 1 );
            bucket_type* offset = table_off(n_large);
            bucket_type* new_s  = offset + (bits_small + 1);
            bucket_type* new_e  = offset + (bits_large + 1);
            std::fill(new_s, new_e, bucket_type());

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

            bucket_type* b0 = table_off(tab);
            bucket_type* b1 = b0 + flag;

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
                        size_type loc = ext::loc(hash, ti);
                        if ( ext::tab(hash, ti) == tab &&
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
        //     auto ntab = std::make_unique<bucket_type[]>(bits_small + 1);
        //     std::vector<std::pair<key_type, mapped_type> > buffer;

        //     migrate_shrnk( n_large, ntab, buffer );

        //     llt[n_large] = std::move(ntab);

        //     finish_shrnk(buffer);

        //     capacity    -= (bits_small+1) * bs;
        //     grow_thresh  = std::ceil((capacity + (bits_large+1)*bs)/alpha);
        //     shrnk_thresh = std::ceil((capacity - (bits_large+1)*bs)/alpha);
        //     if (bits_small == 0 && !n_large) shrnk_thresh = 0;
        // }

        // inline void migrate_shrnk(size_type tab, std::unique_ptr<bucket_type[]>& target,
        //                           std::vector<std::pair<key_type, mapped_type> >& buffer)
        // {
        //     size_type flag = bits_small + 1;

        //     for (size_type i = 0; i < flag; ++i)
        //     {
        //         bucket_type* curr  = &(llt[tab][i]);
        //         bucket_type* curr1 = &(llt[tab][i+flag]);
        //         bucket_type* targ  = &(target[i]);
        //         size_type ind = 0;
        //         for (size_type j = 0; j < bs; ++j)
        //         {
        //             auto e = curr->elements[j];
        //             if (! e.first) break;
        //             auto hash = hasher(e.first);

        //             for (size_type ti = 0; ti < nh; ++ti)
        //             {
        //                 if ( ext::tab(hash, ti)  == tab &&
        //                     (ext::loc(hash, ti) & bits_small) == i)
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
        //                     if ( ext::tab(hash, ti)  == tab &&
        //                          (ext::loc(hash, ti) & bits_small) == i)
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
        //         bla += (base_type::insert(e).second) ? 1: 0;
        //     }
        // }

    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_dysect_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_dysect_inplace<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using key_type         = K;
        using mapped_type      = D;
        using size_type        = size_t;

        static constexpr size_type tl = config_type::tl;
        static constexpr size_type bs = config_type::bs;
        static constexpr size_type nh = config_type::nh;
        static constexpr bool fix_errors = config_type::fix_errors;

        using hasher_type      = hasher<K, HF, ct_log(tl), nh, true, true>;
        using bucket_type      = bucket<K,D,bs>;
    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_dysect_inplace<K,D,HF,Conf> >
    {
    public:
        using table_type = cuckoo_dysect_inplace<K,D,HF,Conf>;
    private:
        using size_type  = typename table_type::size_type;
        using ipointer   = std::pair<const K,D>*;
        static constexpr size_type tl = Conf::tl;
        static constexpr size_type bs = Conf::bs;

    public:
        iterator_incr(const table_type& table_)
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
        const table_type& table;
        ipointer          end_tab;
        size_type         tab;

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

} // namespace dysect
