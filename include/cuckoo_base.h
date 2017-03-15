#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>
#include <limits>

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
    using Config_t       = CuckooConfig<...>;

    using key_type       = ... ;
    using mapped_type    = ... ;

    static constexpr size_t tl = ... ;
    static constexpr size_t bs = ... ;
    static constexpr size_t nh = ... ;

    union Hasher_t       = Hasher<key_type, HashFct, ...>;
    using Bucket_t       = Bucket<key_type, mapped_type, bs>;
};*/

template<class T, bool c = false>
class iterator_incr;
/*
{
    using Table_t = Specialized_t;
    iterator_incr(const Table_t&) { ... }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    pointer next(pointer) { ... }
};
*/


template<class SCuckoo>
class CuckooMultiBase
{
private:
    using  This_t         = CuckooMultiBase<SCuckoo>;
    using  Specialized_t  = typename CuckooTraits<SCuckoo>::Specialized_t;
    using  DisStrat_t     = typename CuckooTraits<SCuckoo>::Config_t::template DisStrat_temp<This_t>;
    using  HistCount_t    = typename CuckooTraits<SCuckoo>::Config_t::HistCount_t;
    using  Bucket_t       = typename CuckooTraits<SCuckoo>::Bucket_t;
    using  Hasher_t       = typename CuckooTraits<SCuckoo>::Hasher_t;
    using  Hashed_t       = typename Hasher_t::Hashed_t;

    friend Specialized_t;
    friend DisStrat_t;

public:
    using key_type        = typename CuckooTraits<SCuckoo>::key_type;
    using mapped_type     = typename CuckooTraits<SCuckoo>::mapped_type;
    using value_type      = std::pair<const key_type, mapped_type>;
    using iterator        = IteratorBase<iterator_incr<Specialized_t> >;
    using const_iterator  = IteratorBase<iterator_incr<Specialized_t>, true>;

private:
    using  value_intern    = std::pair<key_type, mapped_type>;

public:
    CuckooMultiBase(double size_constraint = 1.1,
                    size_t dis_steps = 0, size_t seed = 0);
    ~CuckooMultiBase() = default;
    CuckooMultiBase(const CuckooMultiBase&     ) = delete;
    CuckooMultiBase(      CuckooMultiBase&& rhs);
    CuckooMultiBase& operator=(const CuckooMultiBase&     ) = delete;
    CuckooMultiBase& operator=(      CuckooMultiBase&& rhs)
    {
        n = rhs.n; capacity = rhs.capacity;
        grow_thresh = rhs.grow_thresh; alpha = rhs.alpha;
        return *this;
    }

private:
    // members that should become private at some point ************************
    size_t       n;
    size_t       capacity;
    size_t       grow_thresh;
    double       alpha;
    Hasher_t     hasher;
    DisStrat_t   displacer;
    HistCount_t  hcounter;
    static constexpr size_t bs = CuckooTraits<Specialized_t>::bs;
    static constexpr size_t tl = CuckooTraits<Specialized_t>::tl;
    static constexpr size_t nh = CuckooTraits<Specialized_t>::nh;

public:
    // Basic Hash Table Functionality ******************************************
    iterator                  find  (const key_type& k);
    const_iterator            find  (const key_type& k) const;
    std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
    std::pair<iterator, bool> insert(const value_intern& t);
    size_t                    erase (const key_type& k);

    // Easy use Accessors for std compliance ***********************************
    inline iterator           begin () const; // unimplemented
    inline iterator           end   () const { return make_iterator(nullptr); }
    inline const_iterator     cbegin() const; // unimplemented
    inline const_iterator     cend  () const { return make_citerator(nullptr); }
    mapped_type&              at    (const key_type& k);
    const mapped_type&        at    (const key_type& k) const;
    mapped_type&              operator[](const key_type& k);
    size_t                    count (const key_type& k) const;

private:
    // Easy iterators **********************************************************
    inline iterator           make_iterator (      value_intern* pos) const
    { return iterator      (pos, *static_cast<const Specialized_t*>(this)); }

    inline const_iterator     make_citerator(const value_intern* pos) const
    { return const_iterator(pos, *static_cast<const Specialized_t*>(this)); }

    // implementation specific functions (static polymorph) ********************
    inline void               inc_n() { ++n; }
    inline void               dec_n() { --n; }
    inline void               getBuckets(Hashed_t h, Bucket_t** mem) const
    { return static_cast<const Specialized_t*>(this)->getBuckets(h, mem); }
    inline Bucket_t*          getBucket (Hashed_t h, size_t i) const
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
CuckooMultiBase<SCuckoo>::CuckooMultiBase(double size_constraint,
                                          size_t dis_steps, size_t seed)
    : n(0), capacity(0), grow_thresh(std::numeric_limits<size_t>::max()),
      alpha(size_constraint),
      displacer(*this, dis_steps, seed),
      hcounter(dis_steps)
{ }

template<class SCuckoo>
CuckooMultiBase<SCuckoo>::CuckooMultiBase(CuckooMultiBase&& rhs)
    : n(rhs.n), capacity(rhs.capacity), alpha(rhs.alpha),
      displacer(*this, std::move(rhs.displacer))
{ }



// Implementation of main functionality ****************************************

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::iterator
CuckooMultiBase<SCuckoo>::find(const key_type& k)
{
    auto hash = hasher(k);

    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        value_intern*   tp = tb->findPtr(k);
        if (tp) return make_iterator(tp);
    }
    return end();
}

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::const_iterator
CuckooMultiBase<SCuckoo>::find(const key_type& k) const
{
    auto hash = hasher(k);

    for (size_t i = 0; i < nh; ++i)
    {
        Bucket_t* tb = getBucket(hash, i);
        value_intern*   tp = tb->findPtr(k);
        if (tp) return make_citerator(tp);
    }
    return end();
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const key_type& k, const mapped_type& d)
{
    return insert(std::make_pair(k,d));
}

template<class SCuckoo>
inline std::pair<typename CuckooMultiBase<SCuckoo>::iterator, bool>
CuckooMultiBase<SCuckoo>::insert(const value_intern& t)
{
    if (n > grow_thresh) static_cast<Specialized_t*>(this)->grow();
    auto hash = hasher(t.first);

    std::pair<int,value_intern*> max = std::make_pair(0, nullptr);
    for (size_t i = 0; i < nh; ++i)
    {
        auto temp = getBucket(hash, i)->probePtr(t.first);

        if (temp.first < 0)
            return std::make_pair(make_iterator(temp.second), false);
        max = (max.first > temp.first) ? max : temp;
    }

    if (max.first > 0)
    {
        *max.second = t;
        hcounter.add(0);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(make_iterator(max.second), true);
    }

    int  srch = -1;
    value_intern* pos  = nullptr;
    std::tie(srch, pos) = displacer.insert(t, hash);
    if (srch >=0)
    {
        hcounter.add(srch);
        static_cast<Specialized_t*>(this)->inc_n();
        return std::make_pair(make_iterator(pos), true);
    }

    return std::make_pair(end(), false);
}

template<class SCuckoo>
inline size_t CuckooMultiBase<SCuckoo>::erase(const key_type& k)
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
inline typename CuckooMultiBase<SCuckoo>::mapped_type&
CuckooMultiBase<SCuckoo>::at(const key_type& k)
{
    auto a = static_cast<Specialized_t*>(this)->find(k);
    if (a == end()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SCuckoo>
inline const typename CuckooMultiBase<SCuckoo>::mapped_type&
CuckooMultiBase<SCuckoo>::at(const key_type& k) const
{
    auto a = static_cast<const Specialized_t*>(this)->find(k);
    if (a == cend()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SCuckoo>
inline typename CuckooMultiBase<SCuckoo>::mapped_type&
CuckooMultiBase<SCuckoo>::operator[](const key_type& k)
{
    auto t = static_cast<Specialized_t*>(this)->insert(k, mapped_type());
    return (*t.first).second;
}

template<class SCuckoo>
inline size_t CuckooMultiBase<SCuckoo>::count(const key_type& k) const
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
