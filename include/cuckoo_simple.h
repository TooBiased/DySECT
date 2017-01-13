#pragma once

#include "cuckoo_base.h"

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooSimple : public CuckooTraits<CuckooSimple<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooSimple<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    friend Base_t;

public:
    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;
    static constexpr size_t tl = 1;

private:
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;


    CuckooSimple(size_t cap = 0      , double size_constraint = 1.1,
                      size_t dis_steps = 0, size_t seed = 0)
        : Base_t(std::max(size_t((cap*size_constraint)/bs)*bs, bs), size_constraint,
                 dis_steps, seed),
          beta((size_constraint + 1.)/2.), thresh(cap*beta),
          n_buckets(std::max(size_t((cap*size_constraint)/bs), 1ul)),
          factor(double(n_buckets)/double(1ull<<32)),
          table(new Bucket_t[n_buckets])
    { }

    CuckooSimple(const CuckooSimple&) = delete;
    CuckooSimple& operator=(const CuckooSimple&) = delete;

    CuckooSimple(CuckooSimple&&) = default;
    CuckooSimple& operator=(CuckooSimple&&) = default;

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
                     : std::make_pair(0,nullptr);
    }


    using Base_t::insert;
private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::alpha;
    double beta;
    size_t thresh;
    using Base_t::hasher;

    size_t n_buckets;
    double factor;
    std::unique_ptr<Bucket_t[]> table;
    std::vector<std::pair<Key, Data> > grow_buffer;


    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
        {
            mem[i] = getBucket(h,i);
        }
    }

    inline Bucket_t* getBucket(Hashed_t h, size_t i) const
    {
        size_t l = Ext::loc(h,i) * factor;
        return &(table[l]);
    }

    inline void inc_n()
    {
        if (++n > thresh) grow();
    }

    inline void grow()
    {
        size_t nsize = double(n) * alpha / double(bs);
        nsize        = std::max(nsize, n+1);
        capacity     = nsize*bs;
        double nfactor = double(nsize)/double(1ull << 32);
        size_t nthresh = n * beta;

        auto ntable = std::make_unique<Bucket_t[]>(nsize);
        migrate(ntable, nfactor);
        if (grow_buffer.size()) finalize_grow();

        n_buckets = nsize;
        table     = std::move(ntable);
        thresh    = nthresh;
        factor    = nfactor;
    }

    inline void migrate(std::unique_ptr<Bucket_t[]>& target, double nfactor)
    {
        for (size_t i = 0; i < n_buckets; ++i)
        {
            Bucket_t& curr = table[i];

            for (size_t j = 0; j < bs; ++j)
            {
                auto e = curr.elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);
                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (i == size_t(Ext::loc(hash, i)*factor))
                    {
                        if (! target[Ext::loc(hash, ti) * nfactor].insert(e.first, e.second))
                            grow_buffer.push_back(e);
                        break;
                    }
                }
            }
        }
    }

    inline void finalize_grow()
    {
        size_t temp = n;
        for (auto& e : grow_buffer)
        {
            insert(e);
        }
        n = temp;
        grow_buffer.clear();
    }
};


template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooSimple<K,D,HF,Conf> >
{
public:
    using Specialized_t = CuckooSimple<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using Config_t      = Conf;

    static constexpr size_t tl = 1;
    static constexpr size_t bs = Conf::bs;
    static constexpr size_t nh = Conf::nh;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = Bucket<K,D,bs>;
};
