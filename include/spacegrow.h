#pragma once

#include <functional>
#include <cmath>
#include <memory>
#include <iostream>
#include <vector>
#include <tuple>
#include <random>


#include "bucket.h"
#include "strategies/dstrat_triv.h"

#define MIN_LLS 256

template<class K, class D, class H = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 128, size_t BS = 4>   // ALPHA (SIZE CONSTRAINT COULD BE PARAMETER BUT INTEGRAL TYPE PROBLEM)
class SpaceGrow
{
public:
    using Key       = K;
    using Data      = D;
    using FRet      = std::pair<bool, Data>;

    SpaceGrow(size_t _capacity = 0, double size_constraint = 1.1,
              size_t dis_steps = 0, size_t seed = 0);

    ~SpaceGrow() = default;

    SpaceGrow(const SpaceGrow&) = delete;
    SpaceGrow& operator=(const SpaceGrow&) = delete;

    SpaceGrow(SpaceGrow&& rhs) : nElements(rhs.nElements), curGrowAmount(rhs.curGrowAmount), curGrowTable(rhs.curGrowTable),
                                 capacity(rhs.capacity), alpha(rhs.alpha), bs(BS), tl(TL), displacer(*this, std::move(rhs.displacer))
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
    }

    SpaceGrow& operator=(SpaceGrow&& rhs)
    {
        //if (&rhs == this) return *this;

        //delete this;
        //new (this) SpaceGrow(std::move(rhs));
        //return *this;

        nElements     = rhs.nElements;
        curGrowAmount = rhs.curGrowAmount;
        curGrowTable  = rhs.curGrowTable;
        capacity      = rhs.capacity;
        alpha         = rhs.alpha;

        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = rhs.llb[i];
            llt[i] = std::move(rhs.llt[i]);
        }
        return *this;
    }

    bool insert(Key k, Data d);
    bool insert(std::pair<Key, Data> t);
    FRet find  (Key k);
    bool remove(Key k);

    void clearHist();

private:
    using  This_t     = SpaceGrow<K,D,H,DS,TL,BS>;
    using  Bucket_t   = Bucket<K,D,BS>;
    using  HashFct_t  = H;
    using  DisStrat_t = DS<This_t>;
    friend DisStrat_t;

    static constexpr size_t log(size_t k)
    {
        return (k-1) ? 1+log(k>>1) : 0;
    }

    union HashSplitter {
        std::uint64_t hash;
        struct
        {
            uint64_t tab1 : log(TL);
            uint64_t tab2 : log(TL);
            uint64_t loc1 : 32-log(TL);
            uint64_t loc2 : 32-log(TL);
        };
    };

    static_assert( sizeof(HashSplitter)==8,
                   "HashSplitter must be 64bit!" );

    static_assert( TL == 1ull<<(log(TL)),
                   "TL must be a power of two >0!");

    HashSplitter h(Key k)
    {
        HashSplitter a;
        a.hash = hasher(k);
        return a;
    }

public: //temporary should be removed

    size_t       nElements;
    size_t       curGrowAmount;
    size_t       curGrowTable;
    size_t       capacity;
    double       alpha;
    const size_t bs = BS;
    const size_t tl = TL;
    HashFct_t    hasher;
    DisStrat_t   displacer;

    alignas(64) size_t                      llb[TL];
    alignas(64) std::unique_ptr<Bucket_t[]> llt[TL];  // lower level tables

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        if (i < TL)
            return std::make_pair(llb[i]+1, llt[i].get());
        else
            return std::make_pair(0, nullptr);
    }

    void incElements();
    void grow();
    void migrate(size_t i, std::unique_ptr<Bucket_t[]>& target, size_t tBitmask);

    Bucket_t* getBucket1(HashSplitter h)
    {   return &(llt[h.tab1][(h.loc1 & llb[h.tab1])]);   }
    Bucket_t* getBucket2(HashSplitter h)
    {
        return &(llt[h.tab2][(h.loc2 & llb[h.tab2])]);
    }

};


template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline SpaceGrow<K,D,HF,DS,TL,BS>::SpaceGrow(size_t _capacity, double size_constraint,
                                      size_t dis_steps, size_t seed)
    : nElements(0), alpha(size_constraint), displacer(*this, dis_steps, seed)
{
    std::cout << "this_size:      " << sizeof(This_t)       << std::endl;
    std::cout << "dis_strat_size: " << sizeof(DisStrat_t)   << std::endl;
    std::cout << "hasher_size:    " << sizeof(HashFct_t)    << std::endl;

    double dIni = double(_capacity) * size_constraint / double(TL);
    if (dIni < MIN_LLS)
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = (MIN_LLS/BS) - 1;
            llt[i] = std::make_unique<Bucket_t[]>(MIN_LLS/BS);
        }
        curGrowAmount = MIN_LLS;
        curGrowTable  = 0;
        capacity      = MIN_LLS*TL;

    }
    else
    {
        size_t iIni = MIN_LLS;
        while (dIni > (iIni << 1)) iIni <<= 1;

        size_t gIni = std::floor(double(_capacity) * size_constraint / double(iIni))-TL;

        for (size_t i = gIni; i < TL; ++i)
        {
            llb[i] = (iIni/BS)-1;
            llt[i] = std::make_unique<Bucket_t[]>(iIni/BS);
        }

        curGrowAmount = iIni;
        curGrowTable  = gIni;
        capacity       = (gIni+TL) * iIni;

        iIni         <<= 1;

        for (size_t i = 0; i < gIni; ++i)
        {
            llb[i] = (iIni/BS)-1;
            llt[i] = std::make_unique<Bucket_t[]>(iIni/BS);
        }

    }

    std::cout << "capacity:   " << capacity << std::endl;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(std::pair<Key, Data> t)
{
    auto hash = h(t.first);

    auto p1 = getBucket1(hash)->probe(t.first);//llt[ hash.tab1 ].probe( k, hash.loc1 );
    auto p2 = getBucket2(hash)->probe(t.first);//llt[ hash.tab2 ].probe( k, hash.loc2 );

    if ((p1 < 0) || (p2 < 0)) return false;

    auto r = false;

    // SHOULD CHECK, IF ALREADY INCLUDED
    if (p1 > p2)
    {
        // insert into p1
        r = getBucket1(hash)->insert(t);
    }
    else if (p2 > p1)
    {
        // insert into p2
        r = getBucket2(hash)->insert(t);
    }
    else if (p1)
    {
        // insert into p1
        r = getBucket1(hash)->insert(t);
    }
    else
    {
        // no space => displace stuff
        r = displacer.insert(t, hash);

    }

    if (r) incElements();

    return r;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline typename SpaceGrow<K,D,HF,DS,TL,BS>::FRet SpaceGrow<K,D,HF,DS,TL,BS>::find(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k);
    auto p2 = getBucket2(hash)->find(k);

    if (p1.first) return p1;
    else          return p2;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline bool SpaceGrow<K,D,HF,DS,TL,BS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k);
    auto p2 = getBucket2(hash)->remove(k);

    if (p1) return true;
    else    return p2;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline void SpaceGrow<K,D,HF,DS,TL,BS>::grow()
{
    size_t nsize = (curGrowAmount << 1)/BS; //double the previous size
    auto   ntab  = std::make_unique<Bucket_t[]>(nsize);
    migrate(curGrowTable, ntab, nsize-1);
    llb[curGrowTable] = nsize-1;
    llt[curGrowTable] = std::move(ntab);
    capacity          += curGrowAmount;
    if (++curGrowTable == TL) { curGrowTable = 0; curGrowAmount <<= 1; }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline void SpaceGrow<K,D,HF,DS,TL,BS>::incElements()
{
    ++nElements;
    if (capacity + curGrowAmount < nElements * alpha) grow();
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline void SpaceGrow<K,D,HF,DS,TL,BS>::migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, size_t bitmask)
{
    for (size_t i = 0; i <= llb[tab]; ++i) //Bucket_t* curr = &(llt[i][0]); curr <= &(llt[i][llb[i]]); ++curr)
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

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
inline void SpaceGrow<K,D,HF,DS,TL,BS>::clearHist()
{
    for (size_t i = 0; i < displacer.steps; ++i)
    {
        displacer.hist[i] = 0;
    }
}
