#pragma once

#include "cuckoo_base.h"
#include "strategies/dstrat_triv.h"

union SimpleHashSplitter
{
    uint64_t hash;
    struct
    {
        uint64_t loc1 : 32;
        uint64_t loc2 : 32;
    };
};

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t /*does nothing*/ X = 0, size_t BS = 8>
class SimpleCuckoo : public CuckooBase<K, D, HF, DS, BS,
                                       SimpleCuckoo<K,D,HF,DS,X,BS>,
                                       SimpleHashSplitter>
{
private:
    using This_t       = SimpleCuckoo<K,D,HF,DS,X,BS>;
    using Base_t       = CuckooBase<K,D,HF,DS,BS,This_t,SimpleHashSplitter>;
    friend Base_t;
    using Bucket_t     = typename Base_t::Bucket_t;
    using HashSplitter = typename Base_t::HashSplitter;

public:
    SimpleCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t seed = 0)
        : Base_t(std::max(size_t((cap*size_constraint)/BS)*BS, BS), size_constraint,
                 dis_steps, seed),
          n_buckets(std::max(size_t((cap*size_constraint)/BS), 1ul)),
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

    static constexpr size_t tl = 1;

    using Base_t::n;

private:
    size_t n_buckets;
    std::unique_ptr<Bucket_t[]> table;

    inline Bucket_t* getBucket1(HashSplitter h) const
    { return &(table[h.loc1 % n_buckets]); }

    inline Bucket_t* getBucket2(HashSplitter h) const
    { return &(table[h.loc2 % n_buckets]); }

    inline void inc_n() { ++n; }

};
