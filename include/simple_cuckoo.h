#pragma once

#include "config.h"


template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class SimpleCuckoo : public CuckooTraits<SimpleCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = SimpleCuckoo<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hashed_t       = typename CuckooTraits<This_t>::Hashed_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = CuckooTraits<This_t>::Config_t::bs;
    static constexpr size_t tl = 1;

    SimpleCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t seed = 0)
        : Base_t(std::max(size_t((cap*size_constraint)/bs)*bs, bs), size_constraint,
                 dis_steps, seed),
          n_buckets(std::max(size_t((cap*size_constraint)/bs), 1ul)), factor(double(n_buckets)/double(1ull<<32)),
          table(new Bucket_t[n_buckets])
    { }

    SimpleCuckoo(const SimpleCuckoo&) = delete;
    SimpleCuckoo& operator=(const SimpleCuckoo&) = delete;

    SimpleCuckoo(SimpleCuckoo&&) = default;
    SimpleCuckoo& operator=(SimpleCuckoo&&) = default;

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
                     : std::make_pair(0,nullptr);
    }

private:
    size_t n_buckets;
    double factor;
    std::unique_ptr<Bucket_t[]> table;

    inline Bucket_t* getBucket1(Hashed_t h) const
    { return &(table[h.loc1 * factor]); } //% n_buckets]); }

    inline Bucket_t* getBucket2(Hashed_t h) const
    { return &(table[h.loc2 * factor]); } //% n_buckets]); }
};


template<class K, class D, class HF,
         class Conf>
class CuckooTraits<SimpleCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = SimpleCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using HashFct_t      = HF;
    using Config_t       = Conf;

    static constexpr size_t tl = 1;
    static constexpr size_t bs = Config_t::bs;

    using Bucket_t       = Bucket<K,D,bs>;

    union Hashed_t
    {
        uint64_t hash;
        struct
        {
            uint64_t loc1 : 32;
            uint64_t loc2 : 32;
        };
    };
};
