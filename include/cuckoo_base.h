#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "bucket.h"

// CRTP base class for all cuckoo tables, this encapsulates
// main cuckoo table functionality (insert, find, and remove)
template<class K, class D,
         class HF, template<class> class DS,
         size_t BS, class SCuckoo, class HS>
class CuckooBase
{
private:
    using This_t = CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>;
    using Specialized_t  = SCuckoo;
    using Bucket_t       = Bucket<K,D,BS>;
    using HashFct_t      = HF;
    using DisStrat_t     = DS<This_t>;
    using HashSplitter = HS;
    friend Specialized_t;
    friend DisStrat_t;
    static_assert( sizeof(HashSplitter) == 8,
                   "HashSplitter must be 64bit!" );

public:
    using Key  = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    CuckooBase(size_t cap = 0, double size_constraint = 1.1,
               size_t dis_steps = 0, size_t seed = 0)
        : n(0), capacity(cap), alpha(size_constraint),
          displacer(*this, dis_steps, seed)
    { }

    ~CuckooBase() = default;

    CuckooBase(const CuckooBase&) = delete;
    CuckooBase& operator=(const CuckooBase&) = delete;

    CuckooBase(CuckooBase&& rhs)
        : n(rhs.n), capacity(rhs.capacity), alpha(rhs.alpha),
          displacer(*this, std::move(rhs.displacer))
    { }

    CuckooBase& operator=(CuckooBase&& rhs)
    {
        n = rhs.n;   capacity = rhs.capacity;   alpha = rhs.alpha;
        //displacer = DisStrat_t(*this, std::move(rhs.displacer));
        return *this;
    }

    bool insert(Key k, Data d);
    bool insert(std::pair<Key,Data> t);
    FRet find  (Key k) const;
    bool remove(Key k);

/*** some functions for easier load visualization (not in final product) ******/
    std::pair<size_t, Bucket_t*> getTable(size_t i)
    { return static_cast<Specialized_t*>(this)->getTable(i); }
    void clearHist()
    { for (size_t i = 0; i < displacer.steps; ++i) displacer.hist[i] = 0; }

/*** members that should become private at some point *************************/
    size_t     n;
    size_t     capacity;
    double     alpha;
    HashFct_t  hasher;
    DisStrat_t displacer;
    static constexpr size_t bs = BS;

private:
    inline HashSplitter h(Key k) const
    {
        HashSplitter a;
        a.hash = hasher(k);
        return a;
    }

/*** static polymorph functions ***********************************************/
    inline void inc_n()
    {        static_cast<Specialized_t*>(this)->inc_n(); }

    inline Bucket_t* getBucket1(HashSplitter h) const
    { return static_cast<const Specialized_t*>(this)->getBucket1(h); }

    inline Bucket_t* getBucket2(HashSplitter h) const
    { return static_cast<const Specialized_t*>(this)->getBucket2(h); }
};



/* IMPLEMENTATION *************************************************************/

template<class K, class D,
         class HF, template<class> class DS,
         size_t BS, class SCuckoo, class HS>
inline bool CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class K, class D,
         class HF, template<class> class DS,
         size_t BS, class SCuckoo, class HS>
inline bool CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>::insert(std::pair<Key, Data> t)
{
    auto hash = h(t.first);

    auto p1 = getBucket1(hash)->probe(t.first);
    auto p2 = getBucket2(hash)->probe(t.first);

    if ((p1 < 0) || (p2 < 0)) return false;

    auto r = false;
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

    if (r) inc_n();

    return r;
}

template<class K, class D,
         class HF, template<class> class DS,
         size_t BS, class SCuckoo, class HS>
inline typename CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>::FRet
CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>::find(Key k) const
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k);
    auto p2 = getBucket2(hash)->find(k);

    if (p1.first) return p1;
    else          return p2;
}

template<class K, class D,
         class HF, template<class> class DS,
         size_t BS, class SCuckoo, class HS>
inline bool CuckooBase<K,D,HF,DS,BS,SCuckoo,HS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k);
    auto p2 = getBucket2(hash)->remove(k);

    if (p1) return true;
    else    return p2;
}
