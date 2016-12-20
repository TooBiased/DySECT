#pragma once

#include "config.h"
#include "cuckoo_multi_base.h"

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class SimpleMultiCuckoo : public CuckooTraits<SimpleMultiCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = SimpleMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;
    static constexpr size_t tl = 1;

    SimpleMultiCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t seed = 0)
        : Base_t(std::max(size_t((cap*size_constraint)/bs)*bs, bs), size_constraint,
                 dis_steps, seed),
          n_buckets(std::max(size_t((cap*size_constraint)/bs), 1ul)), factor(double(n_buckets)/double(1ull<<32)),
          table(new Bucket_t[n_buckets])
    { }

    SimpleMultiCuckoo(const SimpleMultiCuckoo&) = delete;
    SimpleMultiCuckoo& operator=(const SimpleMultiCuckoo&) = delete;

    SimpleMultiCuckoo(SimpleMultiCuckoo&&) = default;
    SimpleMultiCuckoo& operator=(SimpleMultiCuckoo&&) = default;

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
                     : std::make_pair(0,nullptr);
    }

private:
    size_t n_buckets;
    double factor;
    std::unique_ptr<Bucket_t[]> table;
    static constexpr size_t u32bitset = (1ull<<32) -1;

    inline void getBuckets(HashSplitter_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
        {
            size_t l = ((h.loc1 + h.loc2*i) & u32bitset) * factor;
            mem[i] = &(table[l]);
        }
    }

    inline Bucket_t* getBucket(HashSplitter_t h, size_t i) const
    {
        size_t l = ((h.loc1 + h.loc2*i) & u32bitset) * factor;
        return &(table[l]);
    }
};


template<class K, class D, class HF,
         class Conf>
class CuckooTraits<SimpleMultiCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = SimpleMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using HashFct_t      = HF;
    using Config_t       = Conf;

    static constexpr size_t tl = 1;
    static constexpr size_t bs = Conf::bs;
    static constexpr size_t nh = Conf::nh;

    using Bucket_t       = Bucket<K,D,bs>;

    union HashSplitter_t
    {
        uint64_t hash;
        struct
        {
            uint64_t loc1 : 32;
            uint64_t loc2 : 32;
        };
    };
};
