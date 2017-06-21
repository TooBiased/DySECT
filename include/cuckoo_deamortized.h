#pragma once

#include "cuckoo_base.h"

template<class K0, class D0, class HF0, class Conf0>
class CuckooIndependentBase;



// *****************************************************************************
// IN PLACE GROWING ************************************************************
// *****************************************************************************

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooDeAmortized : public CuckooTraits<CuckooDeAmortized<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooDeAmortized<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using size_type      = typename CuckooTraits<This_t>::size_type;
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

private:
    using value_intern = std::pair<key_type, mapped_type>;


    static constexpr size_type bs = CuckooTraits<This_t>::bs;
    static constexpr size_type nh = CuckooTraits<This_t>::nh;
    static constexpr size_type tl = 1;

    static constexpr size_type max_size     = 10ull << 30;
    static constexpr size_type min_buckets  = 4096;
    static constexpr size_type grow_step    = 32;

public:
    CuckooDeAmortized(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          beta((size_constraint + 1.)/2.)
    {
        auto temp = reinterpret_cast<Bucket_t*>(operator new (max_size));
        table = std::unique_ptr<Bucket_t[]>(temp);

        // SET CAPACITY, BITMASKS, THRESHOLD
        size_type tcap = size_type(double(cap) * size_constraint / double(bs));
        size_type tbit = min_buckets << 1; while (tbit <= tcap) tbit >>= 1;
        bitmask_large = tbit - 1;
        bitmask_small = bitmask_large >> 1;

        bucket_cutoff = tcap & (~(grow_step-1));
        capacity      = bucket_cutoff*bs;
        grow_thresh   = size_type(double(capacity)/alpha);

        std::fill(table.get(), table.get()+capacity, value_intern());
    }

    CuckooDeAmortized(const CuckooDeAmortized&) = delete;
    CuckooDeAmortized& operator=(const CuckooDeAmortized&) = delete;

    CuckooDeAmortized(CuckooDeAmortized&&) = default;
    CuckooDeAmortized& operator=(CuckooDeAmortized&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    size_type bucket_cutoff;
    size_type bitmask_large;
    size_type bitmask_small;

    std::unique_ptr<value_intern[]> table;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
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

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_type i = 0; i < nh; ++i)
        {
            mem[i] = getBucket(h,i);
        }
    }

    inline Bucket_t* getBucket(Hashed_t h, size_type i) const
    {
        size_type l = Ext::loc(h,i) & bitmask_large;
        if (l > bucket_cutoff) l &= bitmask small;
        return reinterpret_cast<Bucket_t*>(&(table[l*bs]));
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        size_type ncap    = capacity + grow_steps*bs;
        size_type ncutoff = bucket_cutoff + grow_step;
        grow_thresh       = size_type(double(ncap)/alpha);

        std::fill(table.get()+capacity, table.get()+ncap, value_intern());

        migrate(bucket_cutoff, ncutoff);

        capacity      = ncap;
        bucket_cutoff = ncutoff;

        if (ncutoff >= bitmask_large)
        {
            bitmask_small = bitmask_large;
            bitmask_large = (bitmask_large<<1) + 1;
            bucket_cutoff - 0;
        }

    }

    inline void migrate(size_type ocut, size_type ncut)
    {

    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooDeAmortized<K,D,HF,Conf> >
{
public:
    using Specialized_t = CuckooDeAmortized<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Config_t      = Conf;

    using size_type     = size_t;
    using key_type      = K;
    using mapped_type   = D;

    static constexpr size_type tl = 1;
    static constexpr size_type bs = Conf::bs;
    static constexpr size_type nh = Conf::nh;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooDeAmortized<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooDeAmortized<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using pointer   = std::pair<const K,D>*;
    static constexpr size_type bs = Conf::bs;

public:
    iterator_incr(const Table_t& table_)
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
