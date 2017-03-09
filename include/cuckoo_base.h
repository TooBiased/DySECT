#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "bucket.h"
#include "config.h"
#include "hasher.h"
#include "iterator_base.h"
// CRTP base class for all cuckoo tables, this encapsulates
// main cuckoo table functionality (insert, find, and remove)


template<class T>
class CuckooTraits;
/* EXAMPLE IMPLEMENTATION
{
public:
    using Specialized_t  = T;
    using Base_t         = CuckooMultiBase<T>;
    using Key            = ... ;
    using Data           = ... ;
    using Config_t       = CuckooConfig<...>;

    static constexpr size_t tl = ... ;
    static constexpr size_t bs = ... ;
    static constexpr size_t nh = ... ;

    union Hasher_t = Hasher<Key, HashFct, ...>;
    using Bucket_t       = Bucket<Key,Data,bs>;
};*/


template<class SCuckoo>
class CuckooMultiBase
{
private:

    using This_t         = CuckooMultiBase<SCuckoo>;
    using Specialized_t  = typename CuckooTraits<SCuckoo>::Specialized_t;
    using DisStrat_t     = typename CuckooTraits<SCuckoo>::Config_t::template DisStrat_temp<This_t>;
    using HistCount_t    = typename CuckooTraits<SCuckoo>::Config_t::HistCount_t;
    using Hasher_t       = typename CuckooTraits<SCuckoo>::Hasher_t;

    friend Specialized_t;
    friend DisStrat_t;

public:
    using Bucket_t       = typename CuckooTraits<SCuckoo>::Bucket_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;

    using Key      = typename CuckooTraits<SCuckoo>::Key;
    using Data     = typename CuckooTraits<SCuckoo>::Data;
    using Pair_t   = std::pair<Key,Data>;
    using FRet     = std::pair<bool, Data>;
    using Iterator = IteratorBase<This_t>;

    Iterator end() { return Iterator::end(); }


    CuckooMultiBase(size_t cap = 0, double size_constraint = 1.1,
               size_t dis_steps = 0, size_t seed = 0)
        : n(0), capacity(cap), alpha(size_constraint),
          displacer(*this, dis_steps, seed), hcounter(dis_steps)
    { }

    ~CuckooMultiBase() = default;

    CuckooMultiBase(const CuckooMultiBase&) = delete;
    CuckooMultiBase& operator=(const CuckooMultiBase&) = delete;

    CuckooMultiBase(CuckooMultiBase&& rhs)
        : n(rhs.n), capacity(rhs.capacity), alpha(rhs.alpha),
          displacer(*this, std::move(rhs.displacer))
    { }

    CuckooMultiBase& operator=(CuckooMultiBase&& rhs)
    {
        n = rhs.n;   capacity = rhs.capacity;   alpha = rhs.alpha;
        return *this;
    }

    Iterator find(const Key k);
    std::pair<Iterator, bool> insert(const Key k, const Data d);
    std::pair<Iterator, bool> insert(const Pair_t t);

    bool remove(Key k);

    /*** some functions for easier load visualization (not in final product) **/
    std::pair<size_t, Bucket_t*> getTable(size_t i)
    { return static_cast<Specialized_t*>(this)->getTable(i); }
    void clearHist()
    { for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0; }
    inline static void print_init_header(std::ostream& out)
    {
        out.width(6); out << "bsize";
        out.width(6); out << "ntabl";
        out.width(6); out << "nhash";
        out.width(9); out << "f_cap";
        out << std::flush;
    }
    inline void print_init_data(std::ostream& out)
    {
        out.width(6); out << bs;
        out.width(6); out << tl;
        out.width(6); out << nh;
        out.width(9); out << capacity;
        out << std::flush;
    }

    /*** members that should become private at some point *********************/
    size_t       n;
    size_t       capacity;
    double       alpha;
    Hasher_t     hasher;
    DisStrat_t   displacer;
    HistCount_t  hcounter;
    static constexpr size_t bs = CuckooTraits<Specialized_t>::bs;
    static constexpr size_t tl = CuckooTraits<Specialized_t>::tl;
    static constexpr size_t nh = CuckooTraits<Specialized_t>::nh;

private:
    /*** static polymorph functions *******************************************/
    inline void inc_n() { ++n; }
    inline void dec_n() { --n; }

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        return static_cast<const Specialized_t*>(this)->getBuckets(h, mem);
    }

    inline Bucket_t* getBucket(Hashed_t h, size_t i) const
    {
        return static_cast<const Specialized_t*>(this)->getBucket(h, i);
    }
};



/* IMPLEMENTATION *************************************************************/

/* OLD IMPLEMENTATION (PRE ITERATOR)
template<class SCuckoo>
inline bool CuckooMultiBase<SCuckoo>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class SCuckoo>
inline bool CuckooMultiBase<SCuckoo>::insert(std::pair<Key, Data> t)
{
    auto hash = hasher(t.first);
    Bucket_t* ptr[nh];
    getBuckets(hash, ptr);

    int p[nh];
    for (size_t i = 0; i < nh; ++i)
    {
        p[i] = ptr[i]->probe(t.first);
    }

    size_t maxi = 0;
    int maxv = p[0];
    if (p[0] < 0) return false;

    //Seperated from top part to allow pipelining
    for (size_t i = 1; i < nh; ++i)
    {
        if (p[i] < 0) return false;
        if (p[i] > maxv) { maxi = i; maxv = p[i]; }
    }

    auto r = -1;

    // insert into the table with the most space, or trigger displacement
    if (maxv)
    {
        r = (ptr[maxi]->insert(t)) ? 0 : -1;
    }
    else
    {
        std::tie(r, std::ignore) = displacer.insert(t, hash);
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
inline typename CuckooMultiBase<SCuckoo>::FRet
CuckooMultiBase<SCuckoo>::find(Key k) const
{
    auto hash = hasher(k);
    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        FRet      tr = tb->find(k);
        if (tr.first) return tr;
    }
    return std::make_pair(false, 0ull);
}
*/

template<class SCuckoo>
inline bool CuckooMultiBase<SCuckoo>::remove(const Key k)
{
    auto hash = hasher(k);
    Bucket_t* ptr[nh];
    getBuckets(hash, ptr);

    bool p[nh];
    // Make sure this is unrolled and inlined
    for (size_t i = 0; i < nh; ++i)
    {
        p[i] = ptr[i]->remove(k);
    }

    for (size_t i = 0; i < nh; ++i)
    {
        if (p[i])
        {
            static_cast<Specialized_t*>(this)->dec_n();
            return true;
        }
    }
    return false;
}




template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::Iterator CuckooMultiBase<SCuckoo>::find(const Key k)
{
    auto hash = hasher(k);

    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        Pair_t*   tp = tb->findPtr(k);
        if (tp) return Iterator(tp);
    }
    return Iterator::end();
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::Iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const Key k, const Data d)
{
    return insert(std::make_pair(k,d));
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::Iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const Pair_t t)
{
    auto hash = hasher(t.first);

    std::pair<int,Pair_t*> max = std::make_pair(0, nullptr);
    for (size_t i = 0; i < nh; ++i)
    {
        auto temp = getBucket(hash, i)->probePtr(t.first);

        if (temp.first < 0)
            return std::make_pair(Iterator(temp.second), false);
        max = (max.first > temp.first) ? max : temp;
    }

    if (max.first > 0)
    {
        *max.second = t;
        hcounter.add(0);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(Iterator(max.second), true);
    }

    int  srch = -1;
    Pair_t* pos  = nullptr;
    std::tie(srch, pos) = displacer.insert(t, hash);
    if (srch >=0)
    {
        hcounter.add(srch);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(Iterator(pos), true);
    }

    return std::make_pair(end(), false);
}
