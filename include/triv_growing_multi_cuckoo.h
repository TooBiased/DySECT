#pragma once

#include <cmath>
#include "config.h"
#include "cuckoo_multi_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class TrivGrowingMultiCuckoo
    : public CuckooTraits<TrivGrowingMultiCuckoo<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = TrivGrowingMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;

    static constexpr double fac_div = double (1ull << HashSplitter_t::sign_loc);
    static constexpr size_t scnd_lvl_bitmask = 1ull << HashSplitter_t::sign_loc;

    TrivGrowingMultiCuckoo(size_t cap = 0      , double size_constraint = 1.1,
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

    TrivGrowingMultiCuckoo(const TrivGrowingMultiCuckoo&) = delete;
    TrivGrowingMultiCuckoo& operator=(const TrivGrowingMultiCuckoo&) = delete;

    TrivGrowingMultiCuckoo(TrivGrowingMultiCuckoo&& rhs)
        : Base_t(std::move(rhs)), beta(rhs.beta),
          grow_thresh(rhs.grow_thresh), lsize(rhs.lsize), factor(rhs.factor)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    TrivGrowingMultiCuckoo& operator=(TrivGrowingMultiCuckoo&& rhs)
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
    using Base_t::h;

    size_t lsize;
    size_t factor;
    std::unique_ptr<Bucket_t[]> llt[tl];

    inline void getBuckets(HashSplitter_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h, i);
        //return &(llt[h.tab][h.loc1 * factor]);
    } // % lsize]); }

    inline Bucket_t* getBucket (HashSplitter_t h, size_t i) const
    {
        return &(llt[h.tab][((h.loc1+i*h.loc2) & scnd_lvl_bitmask) * factor]);
    } // % lsize]); }

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
                auto hash = h(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (i == ((hash.loc1+ti*hash.loc2) & scnd_lvl_bitmask) * factor)
                    {
                        target[((hash.loc1+ti*hash.loc2) & scnd_lvl_bitmask) * nfactor];
                    }
                }

                // if      (getBucket1(hash) == curr) target[hash.loc1 * factor].insert(e.first, e.second);
                // else if (getBucket2(hash) == curr) target[hash.loc2 * factor].insert(e.first, e.second);
                // else
                // {
                //     std::cout << "something is wrong neither in first, nor second bucket." << std::endl;
                //     exit(64);
                // }
            }
        }
    }
};

template<class K, class D, class HF, class Conf>
class CuckooTraits<TrivGrowingMultiCuckoo<K,D,HF,Conf> >
{
public:
    using Specialized_t  = TrivGrowingMultiCuckoo<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using HashFct_t      = HF;
    using Config_t       = Conf;

    static constexpr size_t bs = Conf::bs;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t nh = Conf::nh;

    using Bucket_t       = Bucket<K,D,bs>;

    union HashSplitter_t
    {
        static constexpr size_t log(size_t k)
        { return (k-1) ? 1+log(k>>1) : 0; }

        static constexpr size_t sign_loc = 32 - log(tl);

        static_assert( tl == 1ull<<log(tl),
                       "TL must be a power of two >0!");

        uint64_t hash;
        struct
        {
            uint64_t unused : log(tl);
            uint64_t tab    : log(tl);
            uint64_t loc1   : 32 - log(tl);
            uint64_t loc2   : 32 - log(tl);
        };
    };
};
