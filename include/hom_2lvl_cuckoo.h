#pragma once

#include <cmath>
#include "cuckoo_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Config = CuckooConfig<> >
class Hom2LvlCuckoo : public CuckooTraits<Hom2LvlCuckoo<K,D,HF,Config> >::Base_t
{
private:
    using This_t         = Hom2LvlCuckoo<K,D,HF,Config>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = CuckooTraits<This_t>::Config_t::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::Config_t::tl;

    Hom2LvlCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                  size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        size_t l2size = std::floor(cap * size_constraint / double(tl*bs));

        for (size_t i = 0; i < tl; ++i)
        {
            llb[i] = l2size;
            llt[i] = std::make_unique<Bucket_t[]>(l2size);
        }

        capacity    = tl * l2size;
    }

    Hom2LvlCuckoo(const Hom2LvlCuckoo&) = delete;
    Hom2LvlCuckoo& operator=(const Hom2LvlCuckoo&) = delete;

    Hom2LvlCuckoo(Hom2LvlCuckoo&& rhs)
        : Base_t(std::move(rhs))
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    Hom2LvlCuckoo& operator=(Hom2LvlCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(llb[i], rhs.llb[i]);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(llb[i], llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::capacity;

private:
    size_t                      llb[tl];
    std::unique_ptr<Bucket_t[]> llt[tl];

    inline Bucket_t* getBucket1(HashSplitter_t h) const
    { return &(llt[h.tab1][h.loc1 % llb[h.tab1]]); }

    inline Bucket_t* getBucket2(HashSplitter_t h) const
    { return &(llt[h.tab2][h.loc2 % llb[h.tab2]]); }
};

template<class K, class D, class HF,
         class Config>
class CuckooTraits<Hom2LvlCuckoo<K,D,HF,Config> >
{
public:
    using Specialized_t  = Hom2LvlCuckoo<K,D,HF,Config>;
    using Base_t         = CuckooBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using HashFct_t      = HF;
    using Config_t       = Config;

    static constexpr size_t tl = Config::tl;
    static constexpr size_t bs = Config::bs;

    using Bucket_t       = Bucket<K,D,bs>;

    union HashSplitter_t
    {
        static constexpr size_t log(size_t k)
        { return (k-1) ? 1+log(k>>1) : 0; }

        static_assert( tl == 1ull<<log(tl),
                       "TL must be a power of two >0!");

        uint64_t hash;
        struct
        {
            uint64_t tab1 : log(tl);
            uint64_t tab2 : log(tl);
            uint64_t loc1 : 32-log(tl);
            uint64_t loc2 : 32-log(tl);
        };
    };
};
