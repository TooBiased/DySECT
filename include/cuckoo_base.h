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

    using Key            = typename CuckooTraits<SCuckoo>::Key;
    using Data           = typename CuckooTraits<SCuckoo>::Data;
    //using key_type       = Key;
    //using mapped_type    = Data;
    //using value_type     = std::pair<const Key, Data>;
    using Pair_t         = std::pair<Key,Data>;
    using iterator       = IteratorBase<This_t>;
    using const_iterator = IteratorBase<This_t, true>;



    CuckooMultiBase(size_t cap = 0, double size_constraint = 1.1,
                    size_t dis_steps = 0, size_t seed = 0);
    ~CuckooMultiBase() = default;
    CuckooMultiBase(const CuckooMultiBase&     ) = delete;
    CuckooMultiBase(      CuckooMultiBase&& rhs);
    CuckooMultiBase& operator=(const CuckooMultiBase&     ) = delete;
    CuckooMultiBase& operator=(      CuckooMultiBase&& rhs)
    {
        n = rhs.n;   capacity = rhs.capacity;   alpha = rhs.alpha;
        return *this;
    }

    // members that should become private at some point ************************
    size_t       n;
    size_t       capacity;
    double       alpha;
    Hasher_t     hasher;
    DisStrat_t   displacer;
    HistCount_t  hcounter;
    static constexpr size_t bs = CuckooTraits<Specialized_t>::bs;
    static constexpr size_t tl = CuckooTraits<Specialized_t>::tl;
    static constexpr size_t nh = CuckooTraits<Specialized_t>::nh;

    // Basic Hash Table Functionality ******************************************
    iterator                  find  (const Key& k);
    const_iterator            find  (const Key& k) const;
    std::pair<iterator, bool> insert(const Key& k, const Data& d);
    std::pair<iterator, bool> insert(const Pair_t& t);
    size_t                    erase (const Key& k);
    //bool                      remove(const Key k) { return erase(k); }

    // Easy use Accessors for std compliance ***********************************
    inline iterator           end   () const { return iterator::end(); }
    inline iterator           begin () const; // unimplemented
    inline const_iterator     cend  () const { return const_iterator::end(); }
    inline const_iterator     cbegin() const; // unimplemented
    Data&                     at    (const Key& k);
    const Data&               at    (const Key& k) const;
    Data&                     operator[](const Key& k);
    size_t                    count (const Key& k) const;

private:
    /*** static polymorph functions *******************************************/
    inline void               inc_n() { ++n; }
    inline void               dec_n() { --n; }
    inline void               getBuckets(Hashed_t h, Bucket_t** mem) const
    { return static_cast<const Specialized_t*>(this)->getBuckets(h, mem); }
    inline Bucket_t*          getBucket(Hashed_t h, size_t i) const
    { return static_cast<const Specialized_t*>(this)->getBucket(h, i); }

public:
    // auxiliary functions for testing *****************************************
    std::pair<size_t, Bucket_t*> getTable(size_t i);
    void                      clearHist();
    void                      print_init_data(std::ostream& out);
    static void               print_init_header(std::ostream& out)
    {
        out.width(6); out << "bsize";
        out.width(6); out << "ntabl";
        out.width(6); out << "nhash";
        out.width(9); out << "f_cap";
        out << std::flush;
    }
};



// Constructors and Appointments ***********************************************

template<class SCuckoo>
CuckooMultiBase<SCuckoo>::CuckooMultiBase(size_t cap, double size_constraint,
                                 size_t dis_steps, size_t seed)
        : n(0), capacity(cap), alpha(size_constraint),
          displacer(*this, dis_steps, seed), hcounter(dis_steps)
{ }

template<class SCuckoo>
CuckooMultiBase<SCuckoo>::CuckooMultiBase(CuckooMultiBase&& rhs)
    : n(rhs.n), capacity(rhs.capacity), alpha(rhs.alpha),
      displacer(*this, std::move(rhs.displacer))
{ }



// Implementation of main functionality ****************************************

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::iterator
CuckooMultiBase<SCuckoo>::find(const Key& k)
{
    auto hash = hasher(k);

    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        Pair_t*   tp = tb->findPtr(k);
        if (tp) return iterator(tp);
    }
    return iterator::end();
}

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::const_iterator
CuckooMultiBase<SCuckoo>::find(const Key& k) const
{
    auto hash = hasher(k);

    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        Pair_t*   tp = tb->findPtr(k);
        if (tp) return const_iterator(tp);
    }
    return const_iterator::end();
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const Key& k, const Data& d)
{
    return insert(std::make_pair(k,d));
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const Pair_t& t)
{
    auto hash = hasher(t.first);

    std::pair<int,Pair_t*> max = std::make_pair(0, nullptr);
    for (size_t i = 0; i < nh; ++i)
    {
        auto temp = getBucket(hash, i)->probePtr(t.first);

        if (temp.first < 0)
            return std::make_pair(iterator(temp.second), false);
        max = (max.first > temp.first) ? max : temp;
    }

    if (max.first > 0)
    {
        *max.second = t;
        hcounter.add(0);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(iterator(max.second), true);
    }

    int  srch = -1;
    Pair_t* pos  = nullptr;
    std::tie(srch, pos) = displacer.insert(t, hash);
    if (srch >=0)
    {
        hcounter.add(srch);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(iterator(pos), true);
    }

    return std::make_pair(end(), false);
}

template<class SCuckoo>
inline size_t CuckooMultiBase<SCuckoo>::erase(const Key& k)
{
    auto hash = hasher(k);
    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        if (tb->remove(k))
        {
            static_cast<Specialized_t*>(this)->dec_n();
            return 1;
        }
    }
    return 0;
}



// Accessor Implementations ****************************************************

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::Data&
CuckooMultiBase<SCuckoo>::at(const Key& k)
{
    auto a = static_cast<Specialized_t*>(this)->find(k);
    if (a == end()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SCuckoo>
inline const typename CuckooMultiBase<SCuckoo>::Data&
CuckooMultiBase<SCuckoo>::at(const Key& k) const
{
    auto a = static_cast<const Specialized_t*>(this)->find(k);
    if (a == cend()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::Data&
CuckooMultiBase<SCuckoo>::operator[](const Key& k)
{
    auto t = static_cast<Specialized_t*>(this)->insert(k, Data());
    return (*t.first).second;
}

template<class SCuckoo>
inline size_t CuckooMultiBase<SCuckoo>::count(const Key& k) const
{
    return (static_cast<const Specialized_t*>(this)->find(k) != cend()) ? 1 : 0;
}



// Print Parameter Functions ***************************************************

/*template<class SCuckoo>
inline static void CuckooMultiBase<SCuckoo>::print_init_header(std::ostream& out)
{
    out.width(6); out << "bsize";
    out.width(6); out << "ntabl";
    out.width(6); out << "nhash";
    out.width(9); out << "f_cap";
    out << std::flush;
}*/

template<class SCuckoo>
inline void CuckooMultiBase<SCuckoo>::print_init_data(std::ostream& out)
{
    out.width(6); out << bs;
    out.width(6); out << tl;
    out.width(6); out << nh;
    out.width(9); out << capacity;
    out << std::flush;
}

template<class SCuckoo>
inline std::pair<size_t, typename CuckooMultiBase<SCuckoo>::Bucket_t*>
CuckooMultiBase<SCuckoo>::getTable(size_t i)
{
    return static_cast<Specialized_t*>(this)->getTable(i);
}

template<class SCuckoo>
inline void CuckooMultiBase<SCuckoo>::clearHist()
{
    for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0;
}
