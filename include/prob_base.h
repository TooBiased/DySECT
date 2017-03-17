#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "config.h"
#include "bucket.h"
#include "iterator_base.h"

template <class T>
class ProbTraits;
template <class T>
class iterator_incr;

template<class SpProb>
class ProbBase
{
private:
    using This_t         = ProbBase<SpProb>;
    using Specialized_t  = typename ProbTraits<SpProb>::Specialized_t;
    using HashFct_t      = typename ProbTraits<SpProb>::HashFct_t;

    friend Specialized_t;
    friend iterator_incr<This_t>;

public:
    using key_type       = typename ProbTraits<SpProb>::key_type;
    using mapped_type    = typename ProbTraits<SpProb>::mapped_type;
    using iterator       = IteratorBase<iterator_incr<This_t> >;
    using const_iterator = IteratorBase<iterator_incr<This_t>, true>;

private:
    using value_intern         = std::pair<key_type,mapped_type>;

public:
    ProbBase(size_t cap, double alpha)
        : alpha(alpha), beta((alpha+1.)/2.), n(0),
          capacity((cap) ? cap*alpha : 2048*alpha),
          thresh  ((cap) ? cap*beta  : 2048*beta)
    {
        if (capacity) table = std::make_unique<value_intern[]>(capacity);
    }

    ~ProbBase() = default;

    ProbBase(const ProbBase&) = delete;
    ProbBase& operator=(const ProbBase&) = delete;

    ProbBase(ProbBase&& rhs) = default;
    ProbBase& operator=(ProbBase&& rhs) = default;

private:
    /*** members that should become private at some point *********************/
    double     alpha;
    double     beta;
    size_t     n;
    size_t     capacity;
    size_t     thresh;
    HashFct_t  hasher;

    std::unique_ptr<value_intern[]> table;

public:
    // Basic Hash Table Functionality ******************************************
    iterator                  find  (const key_type& k);
    const_iterator            find  (const key_type& k) const;
    std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
    std::pair<iterator, bool> insert(const value_intern& t);
    size_t                    erase (const key_type& k);

    // Easy use Accessors for std compliance ***********************************
    inline iterator           begin ()
    {
        auto temp = make_iterator(&table[0]);
        if (!temp->first) temp++;
        return temp;
    }
    inline const_iterator     begin()  const
    {
        return static_cast<Specialized_t*>(this)->cbegin();
    }
    inline const_iterator     cbegin() const
    {
        auto temp = make_citerator(&table[0]);
        if (!temp->first) temp++;
        return temp;
    }
    inline iterator           end   () { return make_iterator(nullptr);  }
    inline const_iterator     end   () const { return static_cast<Specialized_t*>(this)->cend(); }
    inline const_iterator     cend  () const { return make_citerator(nullptr); }

    mapped_type&              at    (const key_type& k);
    const mapped_type&        at    (const key_type& k) const;
    mapped_type&              operator[](const key_type& k);
    size_t                    count (const key_type& k) const;

private:
    // Easy iterators **********************************************************
    inline iterator make_iterator(value_intern* pair) const
    { return iterator(pair, *this); }
    inline const_iterator make_citerator(value_intern* pair) const
    { return const_iterator(pair, *this); }

    // implementation specific functions (static polymorph) ********************
    inline size_t h(key_type k) const
    { return static_cast<const Specialized_t*>(this)->index(hasher(k)); }

    inline void inc_n()
    { if (++n > thresh) static_cast<Specialized_t*>(this)->grow(); }
    inline void dec_n() { --n; }

    // Private helper function *************************************************
    void propagate_remove(size_t origin);

public:
    inline static void print_init_header(std::ostream& out)
    { out.width(9); out << "f_cap"  << " " <<  std::flush;}
    inline void print_init_data  (std::ostream& out)
    { out.width(9); out << capacity << " " << std::flush;}
};



// Implementation of main functionality ****************************************

template<class SpProb>
inline typename ProbBase<SpProb>::iterator
ProbBase<SpProb>::find(const key_type& k)
{
    auto ind = h(k);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<Specialized_t*>(this)->mod(i);
        auto temp = table[ti];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            return make_iterator(&table[ti]);
        }
    }
    return end();
}

template<class SpProb>
inline typename ProbBase<SpProb>::const_iterator
ProbBase<SpProb>::find(const key_type& k) const
{
    auto ind = h(k);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const Specialized_t*>(this)->mod(i);
        auto temp = table[ti];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            return make_citerator(&table[ti]);
        }
    }
    return cend();
}

template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::iterator, bool>
ProbBase<SpProb>::insert(const key_type& k, const mapped_type& d)
{
    return insert(std::make_pair(k,d));
}

template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::iterator, bool>
ProbBase<SpProb>::insert(const value_intern& t)
{
    auto ind = h(t.first);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<Specialized_t*>(this)->mod(i);
        auto temp = table[ti];
        if ( temp.first == t.first)
        {
            return std::make_pair(make_iterator(&table[ti]), false);
        }
        if ( temp.first == 0 )
        {
            table[ti] = t;
            //hcounter.add(i - ind);
            static_cast<SpProb*>(this)->inc_n();
            return std::make_pair(make_iterator(&table[ti]), true);
        }
    }
    return std::make_pair(end(), false);
}

template<class SpProb>
inline size_t ProbBase<SpProb>::erase(const key_type& k)
{
    auto ind = h(k);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<Specialized_t*>(this)->mod(i);
        auto temp = table[ti];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            dec_n();
            table[ti] = std::make_pair(0,0);
            propagate_remove(ti);
            return 1;
        }
    }
    return 0;
}



// Accessor Implementations ****************************************************

template<class SpProb>
inline typename ProbBase<SpProb>::mapped_type&
ProbBase<SpProb>::at(const key_type& k)
{
    auto a = static_cast<Specialized_t*>(this)->find(k);
    if (a == end()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SpProb>
inline const typename ProbBase<SpProb>::mapped_type&
ProbBase<SpProb>::at(const key_type& k) const
{
    auto a = static_cast<const Specialized_t*>(this)->find(k);
    if (a == cend()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SpProb>
inline typename ProbBase<SpProb>::mapped_type&
ProbBase<SpProb>::operator[](const key_type& k)
{
    auto t = static_cast<Specialized_t*>(this)->insert(k, mapped_type());
    return (*t.first).second;
}

template<class SpProb>
inline size_t ProbBase<SpProb>::count(const key_type& k) const
{
    return (static_cast<const Specialized_t*>(this)->find(k) != cend()) ? 1 : 0;
}



// Private Help Function *******************************************************

template<class SpProb>
inline void ProbBase<SpProb>::propagate_remove(const size_t origin)
{
    size_t tempn = n;
    n = 0;
    for (size_t i = origin+1; ; ++i)
    {
        size_t ti = static_cast<const SpProb*>(this)->mod(i);
        auto temp = table[ti];

        if (temp.first == 0) break;

        table[ti] = std::make_pair(0,0);
        insert(temp);
    }
    n = tempn;
}



// Iterator increment **********************************************************

template<class Specialized>
class iterator_incr<ProbBase<Specialized> >
{
public:
    using Table_t    = ProbBase<Specialized>;

private:
    using key_type    = typename Table_t::key_type;
    using mapped_type = typename Table_t::mapped_type;
    using ipointer    = std::pair<const key_type, mapped_type>*;

public:
    iterator_incr(const Table_t& table_)
        : end_ptr(reinterpret_cast<ipointer>
                  (&table_.table[table_.capacity - 1]))
    { }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    ipointer next(ipointer cur)
    {
        while (cur < end_ptr)
        {
            if ((++cur)->first) return cur;
        }
        return nullptr;
    }

private:
    const ipointer end_ptr;
};
