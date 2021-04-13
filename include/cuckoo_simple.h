#pragma once

/*******************************************************************************
 * include/cuckoo_simple.h
 *
 * Requirements: OverAllocation (only for cuckoo_inplace)
 *
 * cuckoo_simple and cuckoo_inplace are implementations of the common
 * d-ary bucket cuckoo hashing scheme. cuckoo_inplace uses
 * overallocation, to grow inplace without full table reallocation.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.hpp"
#include "utils/fastrange.hpp"

#include "cuckoo_base.h"

namespace dysect
{

    template<class K0, class D0, class HF0, class Conf0>
    class cuckoo_adapter_2lvl;

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_standard : public cuckoo_traits<cuckoo_standard<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type      = cuckoo_standard<K,D,HF,Conf>;
        using base_type      = typename cuckoo_traits<this_type>::base_type;
        using bucket_type    = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type    = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type    = typename hasher_type::hashed_type;
        using ext            = typename hasher_type::extractor_type;

        friend base_type;
        friend iterator_incr<this_type>;

        using indep_var = cuckoo_adapter_2lvl<K,D,HF,Conf>;
        friend indep_var;
    public:
        using size_type      = typename cuckoo_traits<this_type>::size_type;
        using key_type       = typename cuckoo_traits<this_type>::key_type;
        using mapped_type    = typename cuckoo_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

    private:
        using value_intern = std::pair<key_type, mapped_type>;

    public:
        cuckoo_standard(size_type cap = 0        , double size_constraint = 1.1,
                        size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed),
              beta((size_constraint + 1.)/2.)
        {
            n_buckets   = size_type(double(cap)*size_constraint)/bs;
            n_buckets   = std::max<size_type>(n_buckets, 256);

            capacity    = n_buckets*bs;
            grow_thresh = beta*std::max<size_type>(256ull, cap);

            //factor      = double(n_buckets)/double(1ull<<32);

            table       = std::make_unique<bucket_type[]>(n_buckets);
        }

        cuckoo_standard(const cuckoo_standard&) = delete;
        cuckoo_standard& operator=(const cuckoo_standard&) = delete;

        cuckoo_standard(cuckoo_standard&&) = default;
        cuckoo_standard& operator=(cuckoo_standard&&) = default;

        ~cuckoo_standard() = default;

    private:
        using base_type::n;
    public:
        using base_type::capacity;
    private:
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        static constexpr size_type bs = cuckoo_traits<this_type>::bs;
        static constexpr size_type nh = cuckoo_traits<this_type>::nh;
        static constexpr size_type tl = 1;
        static constexpr bool fix_errors = cuckoo_traits<this_type>::fix_errors;
        static constexpr size_type min_grow_buckets = 10;

        size_type n_buckets;
        double beta;
        //double factor;

        std::unique_ptr<bucket_type[]> table;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        using base_type::insert;

        std::pair<size_type, bucket_type*> getTable(size_type i)
        {
            return (! i) ? std::make_pair(n_buckets, table.get())
                : std::make_pair(0,nullptr);
        }

        inline iterator begin()
        {
            auto temp = make_iterator(&table[0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

        inline const_iterator cbegin() const
        {
            auto temp = make_citerator(&table[0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

    private:
        // Functions for finding buckets *******************************************

        inline void get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
            {
                mem[i] = get_bucket(h,i);
            }
        }

        inline bucket_type* get_bucket(hashed_type h, size_type i) const
         {
             size_type l = utils_tm::fastrange32(n_buckets, ext::loc(h,i));// * factor;
             return &(table[l]);
         }



        // Size changes (GROWING) **************************************************

        inline void grow()
        {
            size_type nsize = size_type(double(n)*alpha) / bs;
            nsize        = std::max(nsize, n_buckets+min_grow_buckets);
            capacity     = nsize*bs;
            //double nfactor = double(nsize)/double(1ull << 32);
            size_type nthresh = n * beta;

            //std::cout << n << " " << n_buckets << " -> " << nsize << std::endl;

            auto ntable = std::make_unique<bucket_type[]>(nsize);
            std::vector<value_intern> grow_buffer;
            migrate(ntable, nsize, grow_buffer); // factor -> nsize

            n_buckets   = nsize;
            table       = std::move(ntable);
            grow_thresh = nthresh;
            //factor      = nfactor;
            finalize_grow(grow_buffer);
        }

        inline void migrate(std::unique_ptr<bucket_type[]>& target, size_type nsize,
                            std::vector<value_intern>& grow_buffer)
        {
            for (size_type i = 0; i < n_buckets; ++i)
            {
                bucket_type& curr = table[i];

                for (size_type j = 0; j < bs; ++j)
                {
                    auto e = curr.elements[j];
                    if (! e.first) break;
                    auto hash = hasher(e.first);
                    for (size_type ti = 0; ti < nh; ++ti)
                    {
                        if (i == utils_tm::fastrange32(n_buckets, ext::loc(hash, ti)))//*factor))
                        {
                            if (! target[utils_tm::fastrange32(nsize, ext::loc(hash, ti))]
                                         .insert(e.first, e.second))
                            {
                                grow_buffer.push_back(e);
                            }
                            break;
                        }
                    }
                }
            }
        }

        inline void finalize_grow(std::vector<value_intern>& grow_buffer)
        {
            n -= grow_buffer.size(); // n will be fixed by the insertions
            for (auto& e : grow_buffer)
            {
                insert(e);
            }
        }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_standard<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_standard<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using size_type     = size_t;
        using key_type      = K;
        using mapped_type   = D;

        static constexpr size_type tl = 1;
        static constexpr size_type bs = Conf::bs;
        static constexpr size_type nh = Conf::nh;
        static constexpr bool fix_errors = Conf::fix_errors;

        using hasher_type      = hasher<K, HF, 0, nh, true, true>;
        using bucket_type      = bucket<K,D,bs>;
    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_standard<K,D,HF,Conf> >
    {
    public:
        using table_type = cuckoo_standard<K,D,HF,Conf>;
    private:
        using size_type  = typename table_type::size_type;
        using pointer    = std::pair<const K,D>*;
        static constexpr size_type bs = Conf::bs;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<pointer>
                      (&table_.table[table_.n_buckets-1].elements[bs-1]))
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












// *****************************************************************************
// IN PLACE GROWING ************************************************************
// *****************************************************************************

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_standard_inplace : public cuckoo_traits<cuckoo_standard_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type      = cuckoo_standard_inplace<K,D,HF,Conf>;
        using base_type      = typename cuckoo_traits<this_type>::base_type;
        using bucket_type    = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type    = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type    = typename hasher_type::hashed_type;
        using ext            = typename hasher_type::extractor_type;

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
        static constexpr bool fix_errors = cuckoo_traits<this_type>::fix_errors;

        static constexpr size_type max_size     = 10ull << 30;
        static constexpr size_type min_grow_buckets = 10;

    public:
        cuckoo_standard_inplace(size_type cap = 0      , double size_constraint = 1.1,
                      size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed),
              beta((size_constraint + 1.)/2.)
        {
            auto temp = reinterpret_cast<bucket_type*>(operator new (max_size));
            table = std::unique_ptr<bucket_type[]>(temp);

            n_buckets   = size_type(double(cap)*size_constraint)/bs;
            n_buckets   = std::max<size_type>(n_buckets, 256);

            capacity    = n_buckets*bs;
            grow_thresh = beta*std::max<size_type>(256ull, cap);

            //factor      = double(n_buckets)/double(1ull<<32);

            //table       = std::make_unique<bucket_type[]>(n_buckets);
            std::fill(table.get(), table.get()+n_buckets, bucket_type());
        }

        cuckoo_standard_inplace(const cuckoo_standard_inplace&) = delete;
        cuckoo_standard_inplace& operator=(const cuckoo_standard_inplace&) = delete;

        cuckoo_standard_inplace(cuckoo_standard_inplace&&) = default;
        cuckoo_standard_inplace& operator=(cuckoo_standard_inplace&&) = default;

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        size_type n_buckets;
        double    beta;
        //double    factor;

        std::unique_ptr<bucket_type[]> table;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        using base_type::insert;

        std::pair<size_type, bucket_type*> getTable(size_type i)
        {
            return (! i) ? std::make_pair(n_buckets, table.get())
                : std::make_pair(0,nullptr);
        }

        inline iterator begin()
        {
            auto temp = make_iterator(&table[0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

        inline const_iterator cbegin() const
        {
            auto temp = make_citerator(&table[0].elements[0]);
            if (! temp->first) temp++;
            return temp;
        }

    private:
        // Functions for finding buckets *******************************************

        inline void get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
            {
                mem[i] = get_bucket(h,i);
            }
        }

        inline bucket_type* get_bucket(hashed_type h, size_type i) const
        {
            size_type l = utils_tm::fastrange32(n_buckets, ext::loc(h,i));// * factor;
            return &(table[l]);
        }



        // Size changes (GROWING) **************************************************

        inline void grow()
        {
            size_type nsize   = size_type(double(n)*alpha) / bs;
            nsize             = std::max(nsize, n_buckets+min_grow_buckets);
            capacity          = nsize*bs;
            //double    nfactor = double(nsize)/double(1ull << 32);
            size_type nthresh = n * beta;

            std::fill(table.get()+n_buckets, table.get()+nsize, bucket_type());

            std::vector<value_intern> grow_buffer;
            migrate(nsize, grow_buffer);

            n_buckets   = nsize;
            grow_thresh = nthresh;
            //factor      = nfactor;
            finalize_grow(grow_buffer);
        }

        inline void migrate(size_type nsize, std::vector<value_intern>& grow_buffer)
        {
            for (int i = n_buckets; i >= 0; --i)
            {
                bucket_type& curr = table[i];

                for (size_type j = 0; j < bs; ++j)
                {
                    auto e = curr.elements[j];
                    if (! e.first) break;
                    auto hash = hasher(e.first);

                    for (size_type ti = 0; ti < nh; ++ti)
                    {
                        if (i == int(utils_tm::fastrange32(n_buckets, ext::loc(hash, ti))))
                        {
                            int nbucket = utils_tm::fastrange32(nsize, ext::loc(hash,ti));
                            if ((i == nbucket)
                                || (! table[nbucket].insert(e.first, e.second)) )
                            {
                                grow_buffer.push_back(e);
                            }
                            break;
                        }
                    }
                }
                curr = bucket_type();
            }
        }

        inline void finalize_grow(std::vector<value_intern>& grow_buffer)
        {
            n -= grow_buffer.size(); // will be fixed by the following insertions
            for (auto& e : grow_buffer)
            {
                insert(e);
            }
        }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_standard_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_standard_inplace<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using size_type     = size_t;
        using key_type      = K;
        using mapped_type   = D;

        static constexpr size_type tl = 1;
        static constexpr size_type bs = Conf::bs;
        static constexpr size_type nh = Conf::nh;
        static constexpr bool fix_errors = Conf::fix_errors;

        using hasher_type      = hasher<K, HF, 0, nh, true, true>;
        using bucket_type      = bucket<K,D,bs>;

    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_standard_inplace<K,D,HF,Conf> >
    {
    public:
        using table_type = cuckoo_standard_inplace<K,D,HF,Conf>;
    private:
        using size_type  = typename table_type::size_type;
        using pointer    = std::pair<const K,D>*;
        static constexpr size_type bs = Conf::bs;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<pointer>
                      (&table_.table[table_.n_buckets-1].elements[bs-1]))
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



    template<class K, class D, class HF, class Conf>
    class cuckoo_adapter_2lvl
    {
    private:
        using this_type     = cuckoo_adapter_2lvl<K,D,HF,Conf>;
        using subtable_type = cuckoo_standard<K,D,HF,Conf>;
        using hash_function_type  = HF;

        hash_function_type hasher;

    public:
        using key_type       = typename subtable_type::key_type;
        using mapped_type    = typename subtable_type::mapped_type;

        using iterator       = typename subtable_type::iterator;
        using const_iterator = typename subtable_type::const_iterator;

        static constexpr size_t tl   = 256;
        static constexpr size_t bits = tl-1;

    private:
        subtable_type tables[tl];

    public:
        cuckoo_adapter_2lvl(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
        {
            for (size_t i = 0; i < tl; ++i)
            {
                tables[i] = subtable_type(cap/tl, size_constraint);
            }
        }

        inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
        {
            return insert(std::make_pair(k,d));
        }

        inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
        {
            return tables[getInd(t.first)].insert(t);
        }

        inline iterator find(key_type k)
        {
            return tables[getInd(k)].find(k);
        }

        inline size_t erase(key_type k)
        {
            return tables[getInd(k)].erase(k);
        }

        inline iterator begin()              { return tables[0].begin(); }
        inline iterator end()                { return tables[0].end(); }
        inline const_iterator begin()  const { return tables[0].begin(); }
        inline const_iterator end()    const { return tables[0].end(); }
        inline const_iterator cbegin() const { return tables[0].cbegin(); }
        inline const_iterator cend()   const { return tables[0].cend(); }

    private:
        inline size_t getInd(key_type k)
        {
            return hasher(k) & bits;
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(10) << "f_cap" << std::flush;
        }

        inline void print_init_data(otm::output_type& out)
        {
            size_t cap = 0;
            for (size_t i = 0; i < tl; ++i)
                cap += tables[i].capacity;

            out << otm::width(10) << cap << std::flush;
        }
    };
}
