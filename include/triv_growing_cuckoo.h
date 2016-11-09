#pragma once

#include <cmath>
#include "cuckoo_base.h"
#include "strategies/dstrat_triv.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8, class HC = no_hist_count>
class TTrivGrowingCuckoo : public CuckooTraits<TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,HC> >::Base_t
{
private:
    using This_t         = TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,HC>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using HashSplitter_t = typename CuckooTraits<This_t>::HashSplitter_t;

    friend Base_t;

public:
    using Key            = typename CuckooTraits<This_t>::Key;
    using Data           = typename CuckooTraits<This_t>::Data;

    TTrivGrowingCuckoo(size_t cap = 0      , double size_constraint = 1.1,
                  size_t dis_steps = 0, size_t seed = 0)
        : Base_t(0, size_constraint, dis_steps, seed), beta(2.0/(1.0+size_constraint))
    {
        size_t l2size = std::max(std::floor(cap * size_constraint / double(TL*BS)), 1);

        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = l2size;
            llt[i] = std::make_unique<Bucket_t[]>(l2size);
        }

        capacity    = TL * l2size;
        grow_thresh = capacity*beta;
    }

    TTrivGrowingCuckoo(const TTrivGrowingCuckoo&) = delete;
    TTrivGrowingCuckoo& operator=(const TTrivGrowingCuckoo&) = delete;

    TTrivGrowingCuckoo(TTrivGrowingCuckoo&& rhs)
        : Base_t(std::move(rhs))
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    TTrivGrowingCuckoo& operator=(TTrivGrowingCuckoo&& rhs)
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
        return (i < TL) ? std::make_pair(llb[i]+1, llt[i].get())
                        : std::make_pair(0,nullptr);
    }

    using Base_t::n;
    using Base_t::capacity;

private:
    using Base_t::alpha;
    double beta;
    using Base_t::h;

    size_t                      llb[TL];
    std::unique_ptr<Bucket_t[]> llt[TL];

    inline Bucket_t* getBucket1(HashSplitter_t h) const
    { return &(llt[h.tab][h.loc1 & llb[h.tab]]); }

    inline Bucket_t* getBucket2(HashSplitter_t h) const
    { return &(llt[h.tab][h.loc2 & llb[h.tab]]); }

    inline void inc_n()
    {
        ++n;
        if (n > grow_thresh) grow();
    }

    inline void grow()
    {
        size_t nsize = grow_amount << 1;
        auto   ntab  = std::make_unique<Bucket_t[]>(nsize);
        migrate(grow_table, ntab, nsize-1);

        llb[grow_table] = nsize-1;
        llt[grow_table] = std::move(ntab);

        capacity += grow_amount * BS;
        if (++grow_table == TL) { grow_table = 0; grow_amount <<= 1; }
        grow_thresh = std::ceil((capacity + grow_amount*BS)/alpha);
    }

    inline void migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, size_t bitmask)
    {
        for (size_t i = 0; i <= llb[tab]; ++i)
        {
            Bucket_t* curr = &(llt[tab][i]);
            for (size_t j = 0; j < BS; ++j)
            {
                auto e    = curr->elements[j];
                if (! e.first) break;
                auto hash = h(e.first);
                if      (getBucket1(hash) == curr) target[bitmask & hash.loc1].insert(e.first, e.second);
                else if (getBucket2(hash) == curr) target[bitmask & hash.loc2].insert(e.first, e.second);
                else
                {
                    std::cout << "something is wrong neither in first, nor second bucket." << std::endl;
                    exit(64);
                }
            }
        }
    }
};

template<class K, class D, class HF,
         template<class> class DS,
         size_t TL, size_t BS, class HC>
class CuckooTraits<TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,HC> >
{
public:
    using Specialized_t  = TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,HC>;
    using Base_t         = CuckooBase<Specialized_t>;
    using Key            = K;
    using Data           = D;
    using Bucket_t       = Bucket<K,D,BS>;
    using HashFct_t      = HF;
    using DisStrat_t     = DS<Base_t>;
    using HistCount_t    = HC

    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;

    union HashSplitter_t
    {
        static constexpr size_t log(size_t k)
        { return (k-1) ? 1+log(k>>1) : 0; }

        static constexpr size_t uneven(size_t k)
        { return k & 1ull; }

        static constexpr size_t loc_size(size_t k)
        { return 32 - ((log(k)+uneven(k)) >> 1);}

        static_assert( TL == 1ull<<log(TL),
                       "TL must be a power of two >0!");

        uint64_t hash;
        struct
        {
            uint64_t tab  : log(TL);
            //uint64_t tab2 : log(TL);
            uint64_t loc1 : loc_size(TL);
            uint64_t loc2 : loc_size(TL);
        };
    };
};

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using TrivGrowingCuckoo = TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,no_hist_count>

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 256, size_t BS = 8>
using TrivGrowingCuckooHist = TTrivGrowingCuckoo<K,D,HF,DS,TL,BS,hist_count>
