#pragma once

#include <cmath>
#include "cuckoo_multi_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class Hom2LvlMultiCuckoo : public CuckooTraits<Hom2LvlMultiCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = Hom2LvlMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    friend Base_t;

public:
    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;

private:
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    Hom2LvlMultiCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                       size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        size_t l2size = std::floor(double(cap) * size_constraint / double(tl*bs));

        for (size_t i = 0; i < tl; ++i)
        {
            llt[i] = std::make_unique<Bucket_t[]>(l2size);
        }

        capacity    = tl * l2size * bs;
        factor      = double(l2size) / double(1ull << 32 - ct_log(tl));
    }

    Hom2LvlMultiCuckoo(const Hom2LvlMultiCuckoo&) = delete;
    Hom2LvlMultiCuckoo& operator=(const Hom2LvlMultiCuckoo&) = delete;

    Hom2LvlMultiCuckoo(Hom2LvlMultiCuckoo&& rhs)
        : Base_t(std::move(rhs))
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    Hom2LvlMultiCuckoo& operator=(Hom2LvlMultiCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(factor, rhs.factor);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(llb[i], llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::capacity;

private:
    double                      factor;
    std::unique_ptr<Bucket_t[]> llt[tl];

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket(Hashed_t h, size_t i) const
    {
        size_t tab = Ext::tab(h,i);
        return &(llt[tab][Ext::loc(h,i) * factor);
    }
};

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<Hom2LvlMultiCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = Hom2LvlMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Config_t       = Conf;

    static constexpr size_t tl = Config_t::tl;
    static constexpr size_t bs = Config_t::bs;
    static constexpr size_t nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;
};