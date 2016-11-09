#pragma once

#include <cmath>
#include "cuckoo_base.h"
#include "strategies/dstrat_triv.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8, class HC = no_hist_count>
class THom2LvlCuckoo : public CuckooTraits<THom2LvlCuckoo<K,D,HF,DS,TL,BS,HC> >::Base_t
{
private:
    using This_t         = THom2LvlCuckoo<K,D,HF,DS,TL,BS,HC>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    THom2LvlCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                  size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed)
    {
        size_t l2size = std::floor(cap * size_constraint / double(TL*BS));

        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = l2size;
            llt[i] = std::make_unique<Bucket_t[]>(l2size);
        }

        capacity    = TL * l2size;
    }

    THom2LvlCuckoo(const THom2LvlCuckoo&) = delete;
    THom2LvlCuckoo& operator=(const THom2LvlCuckoo&) = delete;

    THom2LvlCuckoo(THom2LvlCuckoo&& rhs)
        : Base_t(std::move(rhs))
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    THom2LvlCuckoo& operator=(THom2LvlCuckoo&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        for (size_t i = 0; i < TL; ++i)
        {
            std::swap(llb[i], rhs.llb[i]);
            std::swap(llt[i], rhs.llt[i]);
        }
        return *this;
    }

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < TL) ? std::make_pair(llb[i], llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::capacity;

private:
    size_t                      llb[TL];
    std::unique_ptr<Bucket_t[]> llt[TL];

    inline Bucket_t* getBucket1(HashSplitter_t h) const
    { return &(llt[h.tab1][h.loc1 % llb[h.tab1]]); }

    inline Bucket_t* getBucket2(HashSplitter_t h) const
    { return &(llt[h.tab2][h.loc2 % llb[h.tab2]]); }
};

template<class K, class D, class HF,
         template<class> class DS,
         size_t TL, size_t BS, class HC>
class CuckooTraits<THom2LvlCuckoo<K,D,HF,DS,TL,BS,HC> >
{
public:
    using Specialized_t  = THom2LvlCuckoo<K,D,HF,DS,TL,BS,HC>;
    using Base_t         = CuckooBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Bucket_t       = Bucket<K,D,BS>;
    using HashFct_t      = HF;
    using DisStrat_t     = DS<Base_t>;
    using HistCount_t    = HC;

    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;

    union HashSplitter_t
    {
        static constexpr size_t log(size_t k)
        { return (k-1) ? 1+log(k>>1) : 0; }

        static_assert( TL == 1ull<<log(TL),
                       "TL must be a power of two >0!");

        uint64_t hash;
        struct
        {
            uint64_t tab1 : log(TL);
            uint64_t tab2 : log(TL);
            uint64_t loc1 : 32-log(TL);
            uint64_t loc2 : 32-log(TL);
        };
    };
};

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using Hom2LvlCuckoo = THom2LvlCuckoo<K,D,HF,DS,TL,BS, no_hist_count>

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using Hom2LvlCuckooHist = THom2LvlCuckoo<K,D,HF,DS,TL,BS, hist_count>;
