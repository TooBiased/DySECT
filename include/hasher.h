#pragma once

/*******************************************************************************
 * include/hasher.h
 *
 * the hasher class is used, to evaluate all hash functions for any
 * given key. From the result it can extract the correct amount of
 * hashed values and split them into appropriate sub parts (subtable
 * number + in-table offset).
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <cstdint>
#include <cstddef>

namespace dysect
{

    static constexpr size_t ct_log(size_t k)
    { return (k-1) ? 1+ct_log(k>>1) : 0; }

    template<size_t t_w, bool dpair>
    struct hash_value_splitter;

    template <class Hashed, size_t tab_width, bool dpair, bool lcomb>
    class  hash_value_extractor;




/* MAIN CLASS *****************************************************************/
    template <class Key, class HFct, size_t TAB_WIDTH, size_t NH,
              bool DPAIR = true, bool LCOMB = true>
    class hasher
    {
    private:
        using this_type = hasher<Key, HFct, TAB_WIDTH, NH, DPAIR, LCOMB>;
        using hash_function_type = HFct;
        static constexpr size_t tab_w   = TAB_WIDTH;
        static constexpr size_t pair_w  = (DPAIR) ? 32 : 64;
        static constexpr size_t loc_w   = pair_w - tab_w;
        static constexpr size_t n_hpair = (LCOMB) ? 2 : NH;
        static constexpr size_t n_hval  = NH;
        static constexpr size_t n_hfct  = (DPAIR) ? (n_hpair>>1)+(1&n_hpair)
            :  n_hpair;
        static_assert( tab_w < pair_w,
                       "TAB_WIDTH has to be smaller than PAIR_WIDTH.");

        using splitter_type = hash_value_splitter<tab_w, DPAIR>;

    public:
        union hashed_type
        {
            uint64_t hash [n_hfct];
            splitter_type split[n_hfct];
        };

        static_assert (sizeof(hashed_type) == n_hfct*8,
                       "Hashhash_value_splitterter Bigger Than Expected");

        using extractor_type = hash_value_extractor<hashed_type, tab_w, DPAIR, LCOMB>;

        /* hasher (itself) ********************************************************/
    private:
        hash_function_type fct[n_hfct];

    public:

        hasher()
        {
            for (size_t i = 0; i < n_hfct; ++i)
            {
                fct[i] = hash_function_type(  2345745572344267838ull +
                                            i*8768656543548765336ull);
            }
        }

        hashed_type operator()(Key k) const
        {
            hashed_type result;
            for (size_t i = 0; i < n_hfct; ++i)
            {
                result.hash[i] = fct[i](k);
            }
            return result;
        }
    };




/* Hash hash_value_splitterter Specializations **********************************************/
    template<size_t t_w>
    struct hash_value_splitter<t_w, true>
    {
        uint64_t tab0 : t_w;
        uint64_t loc0 : 32-t_w;
        uint64_t tab1 : t_w;
        uint64_t loc1 : 32-t_w;
    };

    template<size_t t_w>
    struct hash_value_splitter<t_w, false>
    {
        uint64_t tab  : t_w;
        uint64_t loc  : 64-t_w;
    };

    template<>
    struct hash_value_splitter<0, true>
    {
        static constexpr uint64_t tab0 = 0;
        static constexpr uint64_t tab1 = 0;
        uint64_t loc0 : 32;
        uint64_t loc1 : 32;
    };

    template<>
    struct hash_value_splitter<0, false>
    {
        static constexpr uint64_t tab = 0;
        uint64_t loc  : 64;
    };

/* Extractor Specializations without LCOMB ************************************/
    template <class Hashed, size_t tab_width>
    class hash_value_extractor<Hashed, tab_width, false, false>
    {
    public:
        inline static size_t tab(const Hashed& h, size_t i)
        { return h.split[i].tab; }
        inline static size_t loc(const Hashed& h, size_t i)
        { return h.split[i].loc; }
    };

    template <class Hashed, size_t tab_width>
    class hash_value_extractor<Hashed, tab_width, true, false>
    {
    public:
        inline static size_t tab(const Hashed& h, size_t i)
        { return (i & 1) ? h.split[i>>1].tab1
                : h.split[i>>1].tab0; }
        inline static size_t loc(const Hashed& h, size_t i)
        { return (i & 1) ? h.split[i>>1].loc1
                : h.split[i>>1].loc0; }
    };

/* Extractor Specializations with LCOMB ***************************************/
    template <class Hashed, size_t tab_width>
    class hash_value_extractor<Hashed, tab_width, false, true>
    {
    private:
        static_assert((tab_width != 0) && (tab_width != 64),
                      "Illegal TAB_WIDTH value 0 or 64.");
        static constexpr size_t tab_mask = (1ull<<tab_width     )-1;
        static constexpr size_t loc_mask = (1ull<<(64-tab_width))-1;
    public:
        inline static size_t tab(const Hashed& h, size_t i)
        { return (h.split[0].tab + i*(h.split[1].tab|1)) & tab_mask; }
        inline static size_t loc(const Hashed& h, size_t i)
        { return (h.split[0].loc + i*(h.split[1].loc|1)) & tab_mask; }
    };

    template <class Hashed, size_t tab_width>
    class hash_value_extractor<Hashed, tab_width, true, true>
    {
    private:
        static constexpr size_t tab_mask = (1ull<<tab_width     )-1;
        static constexpr size_t loc_mask = (1ull<<(32-tab_width))-1;
    public:
        inline static size_t tab(const Hashed& h, size_t i)
        { return (h.split[0].tab0 + i*(h.split[0].tab1|1)) & tab_mask; }
        inline static size_t loc(const Hashed& h, size_t i)
        { return (h.split[0].loc0 + i*(h.split[0].loc1|1)) & loc_mask; }
    };

} // namespace dysect
