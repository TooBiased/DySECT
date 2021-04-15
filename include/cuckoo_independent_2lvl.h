#pragma once

/*******************************************************************************
 * include/cuckoo_independent_2lvl.h
 *
 * independent2l is a variant of cuckoo hashing, where the base table
 * is split similar to DySECT, but all buckets connected to one
 * element are placed in the same subtable, therefore, no balancing
 * happens between subtables.  The tables have to grow in small steps
 * (no doubling) in order to hold the given size constraints. During the
 * migration, the size constraint is violated.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <cmath>
#include "utils/default_hash.hpp"
#include "utils/fastrange.hpp"
#include "cuckoo_base.h"


namespace dysect
{

    template<class T>
    class cuckoo_traits;

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_independent_2lvl
        : public cuckoo_traits<cuckoo_independent_2lvl<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type          = cuckoo_independent_2lvl<K,D,HF,Conf>;
        using base_type          = typename cuckoo_traits<this_type>::base_type;
        using bucket_type        = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type        = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type        = typename hasher_type::hashed_type;
        using ext                = typename hasher_type::extractor_type;

        friend base_type;
        friend iterator_incr<this_type>;

    public:
        using size_type          = typename base_type::size_type;
        using key_type           = typename cuckoo_traits<this_type>::key_type;
        using mapped_type        = typename cuckoo_traits<this_type>::mapped_type;
        using iterator           = typename base_type::iterator;
        using const_iterator     = typename base_type::const_iterator;
        using insert_return_type = typename base_type::insert_return_type;

    private:
        using value_intern   = std::pair<key_type, mapped_type>;

    public:
        cuckoo_independent_2lvl(size_type cap = 0      , double size_constraint = 1.1,
                                size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed), beta((1.0+size_constraint)/2.0)
        {
            size_type lsize       = std::floor(cap * size_constraint / double(tl*bs));
            lsize                 = std::max(lsize, 256ul);
            // double    factor      = double(lsize)/fac_div;
            size_type grow_thresh = double(lsize*bs) / beta;

            for (size_type i = 0; i < tl; ++i)
            {
                ll_tab[i]    = std::make_unique<bucket_type[]>(lsize);
                ll_size[i]   = lsize;
                ll_elem[i]   = 0;
                ll_thresh[i] = grow_thresh;
                // ll_factor[i] = factor;
            }

            capacity    = bs * tl * lsize;
        }

        cuckoo_independent_2lvl(const cuckoo_independent_2lvl&) = delete;
        cuckoo_independent_2lvl& operator=(const cuckoo_independent_2lvl&) = delete;

        cuckoo_independent_2lvl(cuckoo_independent_2lvl&& rhs)
            : base_type(std::move(rhs)), beta(rhs.beta)
        {
            for (size_type i = 0; i < tl; ++i)
            {
                ll_tab[i]    = std::move(rhs.ll_tab[i]);
                ll_size[i]   = rhs.ll_size[i];
                ll_elem[i]   = rhs.ll_elem[i];
                ll_thresh[i] = rhs.ll_thresh[i];
                // ll_factor[i] = rhs.ll_factor[i];
            }
        }

        cuckoo_independent_2lvl& operator=(cuckoo_independent_2lvl&& rhs)
        {
            base_type::operator=(std::move(rhs));
            beta            = rhs.beta;

            for (size_type i = 0; i < tl; ++i)
            {
                std::swap(ll_tab[i]   , rhs.ll_tab[i]);
                std::swap(ll_size[i]  , rhs.ll_size[i]);
                std::swap(ll_elem[i]  , rhs.ll_elem[i]);
                std::swap(ll_thresh[i], rhs.ll_thresh[i]);
                // std::swap(ll_factor[i], rhs.ll_factor[i]);
            }
            return *this;
        }

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::alpha;
        using base_type::hasher;

        static constexpr size_type bs = cuckoo_traits<this_type>::bs;
        static constexpr size_type tl = cuckoo_traits<this_type>::tl;
        static constexpr size_type nh = cuckoo_traits<this_type>::nh;
        //static constexpr double    fac_div = double (1ull << (32 - ct_log(tl)));
        static constexpr uint range_shift = 32-ct_log(tl);

        double    beta;
        size_type ll_size  [tl];
        size_type ll_elem  [tl];
        size_type ll_thresh[tl];
        // double    ll_factor[tl];
        std::unique_ptr<bucket_type[]> ll_tab[tl];

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        // Specialized Funcitions (to keep per table counts) ***********************

        inline insert_return_type insert(const key_type k, const mapped_type d)
        {
            return insert(std::make_pair(k,d));
        }

        inline insert_return_type insert(const std::pair<key_type, mapped_type> t)
        {
            auto hash = hasher(t.first);
            size_type ttl = ext::tab(hash, 0);

            auto result = base_type::insert(t);
            if (result.second)
            {
                auto currsize = ++ll_elem[ttl];
                if (currsize > ll_thresh[ttl]) grow_tab(ttl);
            }
            return result;
        }

        size_type erase(const key_type k)
        {
            auto hash     = hasher(k);
            size_type ttl    = ext::tab(hash, 0);
            size_type nk     = base_type::erase(k);
            ll_elem[ttl] -= nk;
            return nk;
        }



        std::pair<size_type, bucket_type*> getTable(size_type i)
        {
            return (i < tl) ? std::make_pair(ll_size[i], ll_tab[i].get())
                : std::make_pair(0,nullptr);
        }

        iterator begin()
        {
            auto temp = make_iterator(&ll_tab[0][0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

        const_iterator cbegin() const
        {
            auto temp = make_citerator(&ll_tab[0][0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

    private:
        // Functions for finding buckets *******************************************

        inline size_type fastrange(size_type cap, size_type h) const
        {
            return (cap*h) >> range_shift;
        }

        inline void get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
                mem[i] = get_bucket(h, i);
        }

        inline bucket_type* get_bucket (hashed_type h, size_type i) const
        {
            size_type tab = ext::tab(h,0);
            // return &(ll_tab[tab][ext::loc(h,i)*ll_factor[tab]]);
            return &(ll_tab[tab][fastrange(ll_size[tab], ext::loc(h,i))]);
        }



        // Size changes (GROWING) **************************************************
        void grow()
        {
            // necessary for base class not implemented since growing is triggered
            // individually for each subtable in this specialization
        }

        inline void grow_tab(size_type tab)
        {
            //otm::out() << "growing in table" << tab << std::endl;
            //TODO make this prettier
            size_type nsize  = std::floor(double(ll_elem[tab]) * alpha / double(bs));
            nsize            = std::max(nsize, ll_size[tab]+1);
            capacity        += (nsize - ll_size[tab])*bs;
            //double nfactor = double(nsize)      / fac_div;
            size_type nthresh = ll_elem[tab] * beta;

            std::vector<value_intern> grow_buffer;

            auto ntable = std::make_unique<bucket_type[]>(nsize);
            migrate(tab, ntable, nsize, grow_buffer);

            ll_tab[tab]    = std::move(ntable);
            ll_size[tab]   = nsize;
            //ll_factor[tab] = nfactor;
            ll_thresh[tab] = nthresh;
            finalize_grow(grow_buffer);
        }

        inline void migrate(size_type tab,
                            std::unique_ptr<bucket_type[]>& target,
                            size_type nsize,
                            std::vector<value_intern>& grow_buffer)
        {
            size_type csize   = ll_size[tab];
            //double cfactor = ll_factor[tab];
            for (size_type i = 0; i < csize; ++i)
            {
                bucket_type& curr = ll_tab[tab][i];

                for (size_type j = 0; j < bs; ++j)
                {
                    auto e    = curr.elements[j];
                    if (! e.first) break;
                    auto hash = hasher(e.first);

                    for (size_type ti = 0; ti < nh; ++ti)
                    {
                        if (i == fastrange(csize, ext::loc(hash, ti)))
                        {
                            if (! target[fastrange(nsize, ext::loc(hash, ti))].insert(e))
                                grow_buffer.push_back(e);
                            break;
                        }
                    }
                }
            }
        }

        inline void finalize_grow(std::vector<value_intern>& grow_buffer)
        {
            for (auto& e : grow_buffer)
            {
                base_type::insert(e);
            }
        }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF, class Conf>
    class cuckoo_traits<cuckoo_independent_2lvl<K,D,HF,Conf> >
    {
    public:
        using specialized_type  = cuckoo_independent_2lvl<K,D,HF,Conf>;
        using base_type         = cuckoo_base<specialized_type>;
        using config_type       = Conf;

        using key_type       = K;
        using mapped_type    = D;

        static constexpr size_t bs = Conf::bs;
        static constexpr size_t tl = Conf::tl;
        static constexpr size_t nh = Conf::nh;
        static_assert(!config_type::fix_errors,
                      "2lvl independent table does not support fix_errors!");
        static constexpr bool fix_errors = false;

        using hasher_type       = hasher<K, HF, ct_log(tl), nh, true, true>;
        using bucket_type       = bucket<K,D,bs>;
    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_independent_2lvl<K,D,HF,Conf> >
    {
    private:
        using ipointer = std::pair<const K,D>*;
        static constexpr size_t tl = Conf::tl;
        static constexpr size_t bs = Conf::bs;

    public:
        using table_type   = cuckoo_independent_2lvl<K,D,HF,Conf>;

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
        size_t            tab;

        ipointer overflow_tab()
        {
            if (++tab >= tl) return nullptr;
            size_t size = table.ll_size[tab] - 1;
            end_tab = reinterpret_cast<ipointer>(&table.ll_tab[tab][size].elements[bs-1]);
            return reinterpret_cast<ipointer>(&table.ll_tab[tab][0].elements[0]);
        }

        void initialize_tab(ipointer ptr)
        {
            for (size_t i = 0; i < tl; ++i)
            {
                ipointer tab_b_ptr = reinterpret_cast<ipointer>
                    (&table.ll_tab[i][0].elements[0]);
                ipointer tab_e_ptr = tab_b_ptr + table.ll_size[i]*bs;

                if (tab_b_ptr <= ptr && ptr < tab_e_ptr)
                {
                    tab = i;
                    end_tab = tab_e_ptr;
                    return;
                }
            }
        }

    };

} // namespace dysect
