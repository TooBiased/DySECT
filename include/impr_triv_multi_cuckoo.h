#pragma once

#include <cmath>
#include "cuckoo_multi_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class IndTableGrowMultiCuckoo
    : public CuckooTraits<IndTableGrowMultiCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = IndTableGrowMultiCuckoo<K,D,HF,Conf>;
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


    static constexpr double fac_div = double (1ull << (32 - ct_log(tl)));

    IndTableGrowMultiCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                           size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed), beta((1.0+size_constraint)/2.0)
    {
        lsize  = std::floor(cap * size_constraint / double(tl*bs));
        lsize  = std::max(lsize, 1ul);
        factor = double(lsize)/fac_div;

        for (size_t i = 0; i < tl; ++i)
        {
            llt[i] = std::make_unique<Bucket_t[]>(lsize);
        }

        capacity    = bs * tl * lsize;
        grow_thresh = capacity / beta;
    }

    IndTableGrowMultiCuckoo(const IndTableGrowMultiCuckoo&) = delete;
    IndTableGrowMultiCuckoo& operator=(const IndTableGrowMultiCuckoo&) = delete;

    IndTableGrowMultiCuckoo(IndTableGrowMultiCuckoo&& rhs)
        : Base_t(std::move(rhs)), beta(rhs.beta),
          grow_thresh(rhs.grow_thresh), lsize(rhs.lsize), factor(rhs.factor)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    IndTableGrowMultiCuckoo& operator=(IndTableGrowMultiCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));
        beta        = rhs.beta;
        grow_thresh = rhs.grow_thresh;
        lsize       = rhs.lsize;
        factor      = rhs.factor;

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(lsize, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::capacity;

private:
    using Base_t::alpha;
    double beta;
    size_t grow_thresh;
    using Base_t::hasher;

    size_t lsize;
    size_t factor;
    std::unique_ptr<Bucket_t[]> llt[tl];

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h, i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_t i) const
    {
        return &(llt[Ext::tab(h,0)][Ext::loc(h,i)*factor]);
    }

    inline void inc_n()
    {
        if ( ++n > grow_thresh ) grow();
    }

    inline void grow()
    {
        size_t nsize   = std::floor(double(n) * alpha / double(tl*bs));
        nsize = std::max(nsize, lsize+1);
        double nfactor = double(nsize)/fac_div;

        for (size_t i = 0; i < tl; ++i)
        {
            auto ntable = std::make_unique<Bucket_t[]>(nsize);
            migrate(i, ntable, nfactor);
            llt[i] = std::move(ntable);
        }

        lsize       = nsize;
        factor      = nfactor;
        capacity    = bs * tl * lsize;
        grow_thresh = capacity / beta;
    }

    inline void migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, double nfactor)
    {
        for (size_t i = 0; i < lsize; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);
            for (size_t j = 0; j < bs; ++j)
            {
                auto e    = curr->elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (i == Ext::loc(hash, ti) * factor)
                    {
                        target[Ext::loc(hash, ti) * nfactor].insert(e.first, e.second);
                    }
                }
            }
        }
    }
};

template<class K, class D, class HF, class Conf>
class CuckooTraits<IndTableGrowMultiCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = IndTableGrowMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Config_t       = Conf;

    static constexpr size_t bs = Conf::bs;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t nh = Conf::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;
};
