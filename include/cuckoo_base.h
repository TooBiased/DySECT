#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "bucket.h"
#include "strategies/dstrat_triv.h"

// CRTP base class for all cuckoo tables, this encapsulates
// main cuckoo table functionality (insert, find, and remove)

class hist_count
{
public:
    hist_count(size_t s) : steps(s), hist(new size_t[s])
    { for (size_t i = 0; i < s; ++i) { hist[i] = 0; } }

    void add(size_t i) { ++hist[i];}

    const size_t steps;
    std::unique_ptr<size_t[]> hist;
};


class no_hist_count
{
public:
    no_hist_count(size_t) { }
    void add(size_t) { }
    static constexpr size_t  steps = 0;
    static constexpr size_t* hist  = nullptr;
};


template<size_t BS = 8, size_t TL = 256,
         template <class> class DisStrat = dstrat_triv,
         class HistCount = no_hist_count>
struct CuckooConfig
{
    static constexpr size_t bs = BS;
    static constexpr size_t tl = TL;

    template <class T>
    using DisStrat_temp = DisStrat<T>;

    using HistCount_t = HistCount;
};

template<class T>
class CuckooTraits;
/* EXAMPLE IMPLEMENTATION
{
public:
    using Specialized_t  = T;
    using Base_t         = CuckooBase<T>;
    using Key            = ... ;
    using Data           = ... ;
    using Bucket_t       = Bucket<Key,Data,bs>;
    using HashFct_t      = ... ;

    using Config_t       = CuckooConfig<...>;

    union HashSplitter
    {
        uint64_t hash;
        struct
        {
            ... //partial hash bits;
        };
    };
};*/


template<class SCuckoo>
class CuckooBase
{
private:

    using This_t         = CuckooBase<SCuckoo>;
    using Specialized_t  = typename CuckooTraits<SCuckoo>::Specialized_t;
    using Bucket_t       = typename CuckooTraits<SCuckoo>::Bucket_t;
    using HashFct_t      = typename CuckooTraits<SCuckoo>::HashFct_t;
    using DisStrat_t     = typename CuckooTraits<SCuckoo>::Config_t::template DisStrat_temp<This_t>;
    using HistCount_t    = typename CuckooTraits<SCuckoo>::Config_t::HistCount_t;
    using HashSplitter_t = typename CuckooTraits<SCuckoo>::HashSplitter_t;

    friend Specialized_t;
    friend DisStrat_t;

    static_assert( sizeof(HashSplitter_t) == 8,
                   "HashSplitter must be 64bit!" );

public:
    using Key      = typename CuckooTraits<SCuckoo>::Key;
    using Data     = typename CuckooTraits<SCuckoo>::Data;
    using FRet     = std::pair<bool, Data>;

    CuckooBase(size_t cap = 0, double size_constraint = 1.1,
               size_t dis_steps = 0, size_t seed = 0)
        : n(0), capacity(cap), alpha(size_constraint),
          displacer(*this, dis_steps, seed), hcounter(dis_steps)
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
        return *this;
    }

    bool insert(Key k, Data d);
    bool insert(std::pair<Key,Data> t);
    FRet find  (Key k) const;
    bool remove(Key k);

    /*** some functions for easier load visualization (not in final product) **/
    std::pair<size_t, Bucket_t*> getTable(size_t i)
    { return static_cast<Specialized_t*>(this)->getTable(i); }
    void clearHist()
    { for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0; }

    /*** members that should become private at some point *********************/
    size_t     n;
    size_t     capacity;
    double     alpha;
    HashFct_t  hasher;
    DisStrat_t displacer;
    HistCount_t  hcounter;
    static constexpr size_t bs = CuckooTraits<Specialized_t>::bs;
    static constexpr size_t tl = CuckooTraits<Specialized_t>::tl;

private:
    inline HashSplitter_t h(Key k) const
    {
        HashSplitter_t a;
        a.hash = hasher(k);
        return a;
    }

    /*** static polymorph functions *******************************************/
    inline void inc_n() { ++n; }
    inline void dec_n() { --n; }

    inline Bucket_t* getBucket1(HashSplitter_t h) const
    { return static_cast<const Specialized_t*>(this)->getBucket1(h); }

    inline Bucket_t* getBucket2(HashSplitter_t h) const
    { return static_cast<const Specialized_t*>(this)->getBucket2(h); }
};



/* IMPLEMENTATION *************************************************************/

template<class SCuckoo>
inline bool CuckooBase<SCuckoo>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class SCuckoo>
inline bool CuckooBase<SCuckoo>::insert(std::pair<Key, Data> t)
{
    auto hash = h(t.first);

    auto p1 = getBucket1(hash)->probe(t.first);
    auto p2 = getBucket2(hash)->probe(t.first);

    if ((p1 < 0) || (p2 < 0)) return false;

    auto r = false;
    if (p1 > p2)
    {
        // insert into p1
        r = (getBucket1(hash)->insert(t)) ? 0 : -1;

    }
    else if (p2 > p1)
    {
        // insert into p2
        r = (getBucket2(hash)->insert(t)) ? 0 : -1;
    }
    else if (p1)
    {
        // insert into p1
        r = (getBucket1(hash)->insert(t)) ? 0 : -1;
    }
    else
    {
        // no space => displace stuff
        r = displacer.insert(t, hash);
    }

    if (r >= 0)
    {
        hcounter.add(r);
        static_cast<Specialized_t*>(this)->inc_n();
        return true;
    }
    else return false;
}

template<class SCuckoo>
inline typename CuckooBase<SCuckoo>::FRet
CuckooBase<SCuckoo>::find(Key k) const
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k);
    auto p2 = getBucket2(hash)->find(k);

    if (p1.first) return p1;
    else          return p2;
}

template<class SCuckoo>
inline bool CuckooBase<SCuckoo>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k);
    auto p2 = getBucket2(hash)->remove(k);

    if (p1 || p2)
    {
        static_cast<Specialized_t*>(this)->dec_n();
        return true;
    }
    else return false;
}
