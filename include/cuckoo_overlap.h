#pragma once

/*******************************************************************************
 * include/cuckoo_overlap.h
 *
 * cuckoo_overlap is an experimental variant of cuckoo hashing, where
 * two different buckets can overlap by some cells. In theory, this
 * achieves higher load factors.  In practice it seems to be slow
 * (ATTENTION NOT OPTIMIZED) additionally it might be
 * counterproductive for DySECT since we want to encourage exchange
 * between different tables.
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
#include "cobucket.h"

namespace dysect
{

    template<class K, class D, class HF = utils_tm::hash_tm::default_hash,
             class Conf = cuckoo_config<> >
    class cuckoo_overlap : public cuckoo_traits<cuckoo_overlap<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type      = cuckoo_overlap<K,D,HF,Conf>;
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

    public:
        cuckoo_overlap(size_type cap = 0        , double size_constraint = 1.1,
                       size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed),
              beta((size_constraint + 1.)/2.)
        {
            n_subbuckets   = size_type(double(cap)*size_constraint)/sbs;
            n_subbuckets   = std::max<size_type>(n_subbuckets, 256);

            capacity    = n_subbuckets*sbs;
            grow_thresh = beta*std::max<size_type>(256ull, cap);

            //factor      = double(n_subbuckets+1-(bs/sbs))/double(1ull<<32);
            acap        = n_subbuckets+1-(bs/sbs);

            table       = std::make_unique<value_intern[]>(n_subbuckets*sbs);
        }

        cuckoo_overlap(const cuckoo_overlap&) = delete;
        cuckoo_overlap& operator=(const cuckoo_overlap&) = delete;

        cuckoo_overlap(cuckoo_overlap&&) = default;
        cuckoo_overlap& operator=(cuckoo_overlap&&) = default;

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        static constexpr size_type sbs = cuckoo_traits<this_type>::sbs;
        static constexpr size_type bs  = cuckoo_traits<this_type>::bs;
        static constexpr size_type nh  = cuckoo_traits<this_type>::nh;
        static constexpr size_type tl  = 1;
        static constexpr bool fix_errors = cuckoo_traits<this_type>::fix_errors;
        static constexpr size_type min_grow_buckets = 10;

        size_type n_subbuckets;
        size_type acap;
        double beta;
//        double factor;

        std::unique_ptr<value_intern[]> table;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        using base_type::insert;

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

        inline void         get_buckets(hashed_type h, bucket_type** mem) const
        {
            for (size_type i = 0; i < nh; ++i)
            {
                mem[i] = get_bucket(h,i);
            }
        }

        inline bucket_type* get_bucket(hashed_type h, size_type i) const
        {
            size_type l = utils_tm::fastrange32(acap, ext::loc(h,i));// * factor;
            return reinterpret_cast<bucket_type*>(&(table[l*sbs]));
        }



        // Size changes (GROWING) **************************************************

        inline void grow()
        {
            //if (!fix_errors && grow_buffer.size()) return; // I think this should not happen
            size_type nsize = size_type(double(n)*alpha) / sbs;
            nsize           = std::max(nsize, n_subbuckets+min_grow_buckets);
            //double nfactor  = double(nsize+1-(bs/sbs))/double(1ull << 32);
            size_type ncap  = nsize+1-(bs/sbs);
            size_type nthresh = n * beta;

            //std::cout << n << " " << n_buckets << " -> " << nsize << std::endl;

            auto ntable     = std::make_unique<value_intern[]>(nsize*sbs);
            std::vector<value_intern> grow_buffer;
            migrate(ntable, ncap, grow_buffer);

            n_subbuckets    = nsize;
            table           = std::move(ntable);
            capacity        = nsize*sbs;
            grow_thresh     = nthresh;
            //factor          = nfactor;
            acap            = ncap;
            finalize_grow(grow_buffer);
        }

        inline void migrate(std::unique_ptr<value_intern[]>& target, size_type ncap,
                            std::vector<value_intern>& grow_buffer)
        {
            for (size_type i = 0; i < capacity; ++i)
            {
                auto e = table[i];
                if (! e.first) continue;
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    size_type off_ti = utils_tm::fastrange32(acap, ext::loc(hash,ti))*sbs;
                    if (i >= off_ti && i < off_ti+bs)
                    {
                        if (! reinterpret_cast<bucket_type*>
                            (&target[utils_tm::fastrange32(ncap, ext::loc(hash, ti))*sbs])
                            ->insert(e))
                        {
                            grow_buffer.push_back(e);
                        }
                    }
                }
            }
        }

        inline void finalize_grow(std::vector<value_intern>& grow_buffer)
        {
            n -= grow_buffer.size(); // the size will be reincremented by insertions
            for (auto& e : grow_buffer)
            {
                insert(e);
            }
        }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_overlap<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_overlap<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using size_type        = size_t;
        using key_type         = K;
        using mapped_type      = D;

        static constexpr size_type tl  = 1;
        static constexpr size_type bs  = Conf::bs;
        static constexpr size_type sbs = Conf::sbs;
        static constexpr size_type nh  = Conf::nh;
        static constexpr bool fix_errors = Conf::fix_errors;

        using hasher_type      = hasher<K, HF, 0, nh, true, true>;
        using bucket_type      = co_bucket<K,D,bs>;

    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_overlap<K,D,HF,Conf> >
    {
    public:
        using table_type   = cuckoo_overlap<K,D,HF,Conf>;
    private:
        using size_type = typename table_type::size_type;
        using pointer   = std::pair<const K,D>*;
        static constexpr size_type bs = Conf::bs;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<pointer>
                      (&table_.table[table_.capacity-1]))
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
    class cuckoo_overlap_inplace : public cuckoo_traits<cuckoo_overlap_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type         = cuckoo_overlap_inplace<K,D,HF,Conf>;
        using base_type         = typename cuckoo_traits<this_type>::base_type;
        using bucket_type       = typename cuckoo_traits<this_type>::bucket_type;
        using hasher_type       = typename cuckoo_traits<this_type>::hasher_type;
        using hashed_type       = typename hasher_type::hashed_type;
        using ext               = typename hasher_type::extractor_type;

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

        static constexpr size_type bs  = cuckoo_traits<this_type>::bs;
        static constexpr size_type sbs = cuckoo_traits<this_type>::sbs;
        static constexpr size_type nh  = cuckoo_traits<this_type>::nh;
        static constexpr size_type tl  = 1;
        static constexpr bool fix_errors = cuckoo_traits<this_type>::fix_errors;
        static constexpr size_type min_grow_buckets = 10;

        static constexpr size_type max_size     = 10ull << 30;

    public:
        cuckoo_overlap_inplace(size_type cap = 0        , double size_constraint = 1.1,
                               size_type dis_steps = 256, size_type seed = 0)
            : base_type(size_constraint, dis_steps, seed),
              beta((size_constraint + 1.)/2.)
        {
            auto temp    = reinterpret_cast<value_intern*>(operator new (max_size));
            table        = std::unique_ptr<value_intern[]>(temp);

            n_subbuckets = size_type(double(cap)*size_constraint)/sbs;
            n_subbuckets = std::max<size_type>(n_subbuckets, 256);

            capacity     = n_subbuckets*sbs;
            grow_thresh  = beta*std::max<size_type>(256ull, cap);

            //factor       = double(n_subbuckets+1-(bs/sbs))/double(1ull<<32);
            acap         = n_subbuckets+1-(bs/sbs);

            std::fill(table.get(), table.get()+capacity, value_intern());
        }

        cuckoo_overlap_inplace(const cuckoo_overlap_inplace&) = delete;
        cuckoo_overlap_inplace& operator=(const cuckoo_overlap_inplace&) = delete;

        cuckoo_overlap_inplace(cuckoo_overlap_inplace&&) = default;
        cuckoo_overlap_inplace& operator=(cuckoo_overlap_inplace&&) = default;

    private:
        using base_type::n;
        using base_type::capacity;
        using base_type::grow_thresh;
        using base_type::alpha;
        using base_type::hasher;

        size_type n_subbuckets;
        size_type acap;
        double    beta;
        // double    factor;

        std::unique_ptr<value_intern[]> table;

        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        using base_type::insert;

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
            size_type l = utils_tm::fastrange32(acap, ext::loc(h,i));// * factor;
            return reinterpret_cast<bucket_type*>(&(table[l*sbs]));
        }



        // Size changes (GROWING) **************************************************

        inline void grow()
        {
            size_type nsize   = size_type(double(n)*alpha) / sbs;
            nsize             = std::max(nsize, n_subbuckets+min_grow_buckets);
            capacity          = nsize*sbs;
            //double    nfactor = double(nsize+1-(bs/sbs))/double(1ull << 32);
            size_type ncap    = nsize+1-(bs/sbs);
            size_type nthresh = n * beta;

            std::fill(table.get()+capacity, table.get()+(nsize*sbs), value_intern());

            std::vector<value_intern> grow_buffer;
            migrate(ncap, grow_buffer);

            n_subbuckets      = nsize;
            grow_thresh       = nthresh;
            // factor            = nfactor;
            acap = ncap;
            finalize_grow(grow_buffer);
        }

        inline void migrate(size_type ncap, std::vector<value_intern>& grow_buffer)
        {
            for (int i = capacity-1; i >= 0; --i)
            {
                auto e = table[i];
                if  (!e.first) continue;
                table[i] = value_intern();
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    size_type off_i = utils_tm::fastrange32(acap, ext::loc(hash, ti))*sbs;
                    if (i >= int(off_i) && i < int(off_i + bs))
                    {
                        size_type nbucket = utils_tm::fastrange32(ncap, ext::loc(hash,ti))*sbs;
                        if (nbucket >= size_type(i) ||
                            ! reinterpret_cast<bucket_type*>(&table[nbucket])->insert(e))
                        {
                            grow_buffer.push_back(e);
                        }
                        break;
                    }
                }
            }
        }

        inline void finalize_grow(std::vector<value_intern>& grow_buffer)
        {
            n -= grow_buffer.size(); // n will be increased by insertions
            for (auto& e : grow_buffer)
            {
                insert(e);
            }
        }
    };



// Traits class defining types *************************************************

    template<class K, class D, class HF,
             class Conf>
    class cuckoo_traits<cuckoo_overlap_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type = cuckoo_overlap_inplace<K,D,HF,Conf>;
        using base_type        = cuckoo_base<specialized_type>;
        using config_type      = Conf;

        using size_type        = size_t;
        using key_type         = K;
        using mapped_type      = D;

        static constexpr size_type tl = 1;
        static constexpr size_type bs = Conf::bs;
        static constexpr size_type nh = Conf::nh;
        static constexpr size_type sbs= Conf::sbs;
        static constexpr bool fix_errors = Conf::fix_errors;

        using hasher_type      = hasher<K, HF, 0, nh, true, true>;
        using bucket_type      = co_bucket<K,D,bs>;

    };



// Iterator increment **********************************************************

    template<class K, class D, class HF, class Conf>
    class iterator_incr<cuckoo_overlap_inplace<K,D,HF,Conf> >
    {
    public:
        using table_type   = cuckoo_overlap_inplace<K,D,HF,Conf>;
    private:
        using size_type = typename table_type::size_type;
        using pointer   = std::pair<const K,D>*;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<pointer>
                      (&table_.table[table_.capacity-1]))
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
