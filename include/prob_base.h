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

template<class SpProb>
class ProbBase
{
private:
    using This_t         = ProbBase<SpProb>;
    using Specialized_t  = typename ProbTraits<SpProb>::Specialized_t;
    friend Specialized_t;
    using HashFct_t      = typename ProbTraits<SpProb>::HashFct_t;

public:
    using Key      = typename ProbTraits<SpProb>::Key;
    using Data     = typename ProbTraits<SpProb>::Data;
    using Pair_t   = std::pair<Key,Data>;
    using iterator = IteratorBase<This_t>;
    using const_iterator = IteratorBase<This_t, true>;
    //static constexpr iterator end() { return iterator::end(); }


    ProbBase(size_t cap, double alpha)
        : alpha(alpha), beta((alpha+1.)/2.), n(0),
          capacity((cap) ? cap*alpha : 2048*alpha),
          thresh  ((cap) ? cap*beta  : 2048*beta)
    {
        if (capacity) table = std::make_unique<Pair_t[]>(capacity);
    }

    ~ProbBase() = default;

    ProbBase(const ProbBase&) = delete;
    ProbBase& operator=(const ProbBase&) = delete;

    ProbBase(ProbBase&& rhs) = default;
    ProbBase& operator=(ProbBase&& rhs) = default;

    // Basic Hash Table Functionality ******************************************
    iterator                  find  (const Key& k);
    const_iterator            find  (const Key& k) const;
    std::pair<iterator, bool> insert(const Key& k, const Data& d);
    std::pair<iterator, bool> insert(const Pair_t& t);
    size_t                    erase (const Key& k);

    // Easy use Accessors for std compliance ***********************************
    inline iterator           end   () const { return iterator::end(); }
    inline iterator           begin () const; // unimplemented
    inline const_iterator     cend  () const { return const_iterator::end(); }
    inline const_iterator     cbegin() const; // unimplemented
    Data&                     at    (const Key& k);
    const Data&               at    (const Key& k) const;
    Data&                     operator[](const Key& k);
    size_t                    count (const Key& k) const;


    inline static void print_init_header(std::ostream& out)
    { out.width(9); out << "f_cap"  << " " <<  std::flush;}
    inline void print_init_data  (std::ostream& out)
    { out.width(9); out << capacity << " " << std::flush;}

    /*** members that should become private at some point *********************/
    double     alpha;
    double     beta;
    size_t     n;
    size_t     capacity;
    size_t     thresh;
    HashFct_t  hasher;

    std::unique_ptr<Pair_t[]> table;

private:
    /*** static polymorph functions *******************************************/
    inline size_t h(Key k) const
    { return static_cast<const Specialized_t*>(this)->index(hasher(k)); }

    void propagate_remove(size_t origin);

    inline void inc_n()
    { if (++n > thresh) static_cast<Specialized_t*>(this)->grow(); }
    inline void dec_n() { --n; }
};



// Implementation of main functionality ****************************************

template<class SpProb>
inline typename ProbBase<SpProb>::iterator
ProbBase<SpProb>::find(const Key& k)
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
            return iterator(&table[ti]);
        }
    }
    return end();
}

template<class SpProb>
inline typename ProbBase<SpProb>::const_iterator
ProbBase<SpProb>::find(const Key& k) const
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
            return const_iterator(&table[ti]);
        }
    }
    return cend();
}

template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::iterator, bool>
ProbBase<SpProb>::insert(const Key& k, const Data& d)
{
    return insert(std::make_pair(k,d));
}

template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::iterator, bool>
ProbBase<SpProb>::insert(const Pair_t& t)
{
    auto ind = h(t.first);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<Specialized_t*>(this)->mod(i);
        auto temp = table[ti];
        if ( temp.first == t.first)
        {
            return std::make_pair(iterator(&table[ti]), false);
        }
        if ( temp.first == 0 )
        {
            table[ti] = t;
            //hcounter.add(i - ind);
            static_cast<SpProb*>(this)->inc_n();
            return std::make_pair(iterator(&table[ti]), true);
        }
    }
    return std::make_pair(end(), false);
}

template<class SpProb>
inline size_t ProbBase<SpProb>::erase(const Key& k)
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
inline typename ProbBase<SpProb>::Data&
ProbBase<SpProb>::at(const Key& k)
{
    auto a = static_cast<Specialized_t*>(this)->find(k);
    if (a == end()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SpProb>
inline const typename ProbBase<SpProb>::Data&
ProbBase<SpProb>::at(const Key& k) const
{
    auto a = static_cast<const Specialized_t*>(this)->find(k);
    if (a == cend()) throw std::out_of_range("cannot find key");
    else return (*a).second;
}

template<class SpProb>
inline typename ProbBase<SpProb>::Data&
ProbBase<SpProb>::operator[](const Key& k)
{
    auto t = static_cast<Specialized_t*>(this)->insert(k, Data());
    return (*t.first).second;
}

template<class SpProb>
inline size_t ProbBase<SpProb>::count(const Key& k) const
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
