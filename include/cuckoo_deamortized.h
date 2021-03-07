#pragma once

/*******************************************************************************
 * include/cuckoo_deamortized.h
 *
 * Requirement: <OverAllocation>
 *
 * cuckoo_deamortized is a dynamically growing cuckoo table variation
 * that uses only one table. The idea is to grow the table in place
 * but do it in small steps. Therefore, grown parts of the table
 * coexist with ungrown ones.
 *
 * 0          x                 2^k     cap = 2^k+x
 * +----------+------------------+----------+ - - - - - - - - - -
 * | grown    | ungrown          | new part |  overallocated ...
 * +----------+------------------+----------+ - - - - - - - - - -
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.hpp"
#include "cuckoo_base.h"

namespace dysect
{

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_deamortized : public cuckoo_traits<cuckoo_deamortized<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type   = cuckoo_deamortized<K,D,HF,Conf>;
        using base_type   = typename cuckoo_traits<this_type>::base_type;
        using bucket_type = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type = typename hasher_type::hashed_type;
        using ext         = typename hasher_type::extractor_type;

        friend base_type;
        friend iterator_incr<this_type>;

    public:
        using size_type      = typename cuckoo_traits<this_type>::size_type;
        using key_type       = typename cuckoo_traits<this_type>::key_type;
        using mapped_type    = typename cuckoo_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

    private:
        using value_intern = std::pair<key_type, mapped_type>;


        static constexpr size_type bs = cuckoo_traits<this_type>::bs;
        static constexpr size_type nh = cuckoo_traits<this_type>::nh;
        static constexpr size_type tl = 1;

        static constexpr size_type max_size     = 10ull << 30;
        static constexpr size_type min_buckets  = 4096;
        static constexpr size_type grow_step    = 32;

    public:
        cuckoo_deamortized(size_type cap = 0      , double size_constraint = 1.1,
                          size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed)
            {
                auto temp = reinterpret_cast<value_intern*>(operator new (max_size));
                table = std::unique_ptr<value_intern[]>(temp);

                // SET CAPACITY, BITMASKS, THRESHOLD
                size_type tcap = size_type(double(cap) * size_constraint / double(bs));
                size_type tbit = min_buckets << 1; while (tbit <= tcap) tbit <<= 1;
                bitmask_large = tbit - 1;
                bitmask_small = bitmask_large >> 1;

                size_type doub = tcap -bitmask_small-1;
                doub /= grow_step;
                doub *= grow_step;

                capacity      = (doub+bitmask_small+1)*bs;
                bucket_cutoff = doub+bitmask_small+1;
                grow_thresh   = size_type(double(capacity+grow_step*bs)/alpha);

                std::fill(table.get(), table.get()+capacity, value_intern());
            }

        cuckoo_deamortized(const cuckoo_deamortized&) = delete;
        cuckoo_deamortized& operator=(const cuckoo_deamortized&) = delete;

        cuckoo_deamortized(cuckoo_deamortized&&) = default;
        cuckoo_deamortized& operator=(cuckoo_deamortized&&) = default;

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        size_type bucket_cutoff;
        size_type bitmask_large;
        size_type bitmask_small;

        std::unique_ptr<value_intern[]> table;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        using base_type::insert;

        std::pair<size_type, bucket_type*> getTable(size_type i)
            {
                return (! i) ? std::make_pair(bucket_cutoff, table.get())
                    : std::make_pair(0,nullptr);
            }

        inline iterator begin()
            {
                auto temp = make_iterator(&table[0]);
                if (! temp->first) temp++;
                return temp;
            }

        inline const_iterator cbegin() const
            {
                auto temp = make_citerator(&table[0]);
                if (! temp->first) temp++;
                return temp;
            }

    private:
        // Functions for finding buckets *******************************************

        inline void getBuckets(hashed_type h, bucket_type** mem) const
            {
                for (size_type i = 0; i < nh; ++i)
                {
                    mem[i] = getBucket(h,i);
                }
            }

        inline bucket_type* getBucket(hashed_type h, size_type i) const
            {
                size_type l = ext::loc(h,i) & bitmask_large;
                if (l >= bucket_cutoff) l &= bitmask_small;
                return reinterpret_cast<bucket_type*>(&(table[l*bs]));
            }



        // Size changes (GROWING) **************************************************

        inline void grow()
            {
                size_type ncap    = capacity + grow_step*bs;
                size_type ncutoff = bucket_cutoff + grow_step;
                grow_thresh       = size_type(double(ncap+grow_step*bs)/alpha);

                std::fill(table.get()+capacity, table.get()+ncap, value_intern());

                migrate(capacity, ncap);

                capacity      = ncap;
                bucket_cutoff = ncutoff;

                if (ncutoff >= bitmask_large)
                {
                    bitmask_small = bitmask_large;
                    bitmask_large = (bitmask_large<<1) + 1;
                }

            }

        inline void migrate(size_type , size_type)
            {
                size_type flag  = bitmask_large  & (~bitmask_small);
                size_type ptr0  = (bucket_cutoff & bitmask_small) * bs;

                for (size_type i = 0; i < grow_step; ++i)
                {
                    bucket_type* bucket0_ptr =
                        reinterpret_cast<bucket_type*>(table.get()+ptr0+i*bs);
                    bucket_type* bucket1_ptr =
                        reinterpret_cast<bucket_type*>(table.get()+capacity+i*bs);

                    size_type c0 = 0;
                    size_type c1 = 0;
                    for (size_type j = 0; j < bs; ++j)
                    {
                        auto curr = bucket0_ptr->elements[j];
                        bucket0_ptr->elements[j] = value_intern();
                        if (!curr.first) break;

                        bucket_type* targets[nh];
                        hashed_type  hash = hasher(curr.first);
                        getBuckets(hash, targets);

                        for (size_type t = 0; t < nh; ++t)
                        {
                            if (targets[t] == bucket0_ptr)
                            {
                                if (ext::loc(hash, t) & flag)
                                {
                                    bucket1_ptr->elements[c1++] = curr;
                                }
                                else
                                {
                                    bucket0_ptr->elements[c0++] = curr;
                                }
                                break;
                            }
                        }
                    }
                }
            }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_deamortized<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_deamortized<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using size_type        = size_t;
        using key_type         = K;
        using mapped_type      = D;

        static constexpr size_type tl = 1;
        static constexpr size_type bs = Conf::bs;
        static constexpr size_type nh = Conf::nh;
        static constexpr bool fix_errors = Conf::fix_errors;

        using hasher_type      = hasher<K, HF, 0, nh, true, true>;
        using bucket_type      = bucket<K,D,bs>;

    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_deamortized<K,D,HF,Conf> >
    {
    public:
        using table_type = cuckoo_deamortized<K,D,HF,Conf>;
    private:
        using size_type  = typename table_type::size_type;
        using pointer    = std::pair<const K,D>*;
        static constexpr size_type bs = Conf::bs;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<pointer>(&table_.table[table_.capacity -1]))
            { }
        iterator_incr(const iterator_incr&) = default;
        iterator_incr& operator=(const iterator_incr&) = default;

        pointer next(pointer cur)
            {
                while (cur < end_ptr)
                {
                    if ((++cur)->first) return cur;
                }
                return nullptr;
            }

    private:
        pointer end_ptr;
    };

} // namespace dysect
