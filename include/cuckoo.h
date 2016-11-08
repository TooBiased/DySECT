#pragma once

#include <memory>

#include "bucket.h"
#include "strategies/dstrat_triv.h"

template<class K, class D, class HF = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t /*does nothing*/X = 0, size_t BS = 4>
class CuckooTable
{
public:
    using Key = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    CuckooTable(size_t capacity  = 0, double size_constraint = 1.1,
                size_t dis_steps = 0, size_t seed = 0);

    ~CuckooTable() = default;

    CuckooTable(const CuckooTable&)          = delete;
    CuckooTable& operator=(const CuckooTable&) = delete;

    CuckooTable(CuckooTable&&) = default;
    CuckooTable& operator=(CuckooTable&&) = default;

    bool insert(Key k, Data d);
    bool insert(std::pair<Key, Data> t);
    FRet find(Key k);
    bool remove(Key k);

    void clearHist();

private:
    using This_t = CuckooTable<K,D,HF,DS,X,BS>;
    using Bucket_t = Bucket<K,D,BS>;
    using HashFct_t= HF;
    using DisStrat_t = DS<This_t>;
    friend DisStrat_t;

    union HashSplitter {
        std::uint64_t hash;
        struct
        {
            uint64_t loc1 : 32;
            uint64_t loc2 : 32;
        };
    };

    static_assert( sizeof(HashSplitter)==8,
                    "HashSplitter must be 64bit!");
    HashSplitter h(Key k)
    {
        HashSplitter a;
        a.hash = hasher(k);
        return a;
    }

public: //temporary should be removed

    size_t       nElements;
    size_t       nBuckets;
    size_t       capacity;
    double       alpha;
    const size_t bs = BS;
    const size_t tl = 1;
    HashFct_t    hasher;
    DisStrat_t   displacer;

    std::unique_ptr<Bucket_t[]> table;

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        if (i == 0)
            return std::make_pair(nBuckets, table.get());
        else
            return std::make_pair(0, nullptr);
    }

    void incElements();
    //void grow();
    //void migrate(std::unique_ptr<Bucket[]>& target, size_t newNBuckets);

    Bucket_t* getBucket1(HashSplitter h)
    {   return &(table[h.loc1 % nBuckets]);   }
    Bucket_t* getBucket2(HashSplitter h)
    {   return &(table[h.loc2 % nBuckets]);   }
};

template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
CuckooTable<K,D,HF,DS,X,BS>::CuckooTable(size_t cap, double size_constraint,
                                    size_t dis_steps, size_t seed)
    : nElements(0), nBuckets(std::max(size_t((cap*size_constraint)/BS), size_t(256))),
      capacity(nBuckets*BS), alpha(size_constraint),
      displacer(*this, dis_steps, seed), table(new Bucket_t[nBuckets])
{
}

template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline bool CuckooTable<K,D,HF,DS,X,BS>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}


template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline bool CuckooTable<K,D,HF,DS,X,BS>::insert(std::pair<Key, Data> t)
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
    else if (nElements < capacity)
    {
        // no space => displace stuff
        r = displacer.insert(t, hash);
    }

    if (r) incElements();

    return r;
}


template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline typename CuckooTable<K,D,HF,DS,X,BS>::FRet CuckooTable<K,D,HF,DS,X,BS>::find(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k);
    auto p2 = getBucket2(hash)->find(k);

    if (p1.first) return p1;
    else          return p2;
}


template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline bool CuckooTable<K,D,HF,DS,X,BS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k);
    auto p2 = getBucket2(hash)->remove(k);

    if (p1) return true;
    else    return p2;
}


template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline void CuckooTable<K,D,HF,DS,X,BS>::clearHist()
{
    for (size_t i = 0; i < displacer.steps; ++i)
    {
        displacer.hist[i] = 0;
    }
}

template<class K, class D, class HF,
         template<class> class DS,
         size_t X, size_t BS>
inline void CuckooTable<K,D,HF,DS,X,BS>::incElements()
{
    ++nElements;
    //if (capacity + curGrowAmount < nElements * alpha) grow();
}
