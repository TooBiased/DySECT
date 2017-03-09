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
    using Iterator = IteratorBase<This_t>;
    static constexpr Iterator end() { return Iterator::end(); }


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

    std::pair<Iterator, bool> insert(Key k, Data d);
    std::pair<Iterator, bool> insert(std::pair<Key,Data> t);
    Iterator find  (Key k) const;
    bool remove(Key k);// { return false; }

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



/* IMPLEMENTATION *************************************************************/
template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::Iterator, bool>
ProbBase<SpProb>::insert(const Key k, const Data d)
{
    return insert(std::make_pair(k,d));
}

template<class SpProb>
inline std::pair<typename ProbBase<SpProb>::Iterator, bool>
ProbBase<SpProb>::insert(const std::pair<Key, Data> t)
{
    auto ind = h(t.first);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const SpProb*>(this)->mod(i);
        auto temp = table[ti];
        if ( temp.first == t.first)
        {
            return std::make_pair(Iterator(&table[ti]), false);
        }
        if ( temp.first == 0 )
        {
            table[ti] = t;
            //hcounter.add(i - ind);
            static_cast<SpProb*>(this)->inc_n();
            return std::make_pair(Iterator(&table[ti]), true);
        }
    }
    return std::make_pair(end(), false);
}

template<class SpProb>
inline typename ProbBase<SpProb>::Iterator
ProbBase<SpProb>::find(const Key k) const
{
    auto ind = h(k);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const SpProb*>(this)->mod(i);
        auto temp = table[ti];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            return Iterator(&table[ti]);
        }
    }
    return end();
}

template<class SpProb>
inline bool ProbBase<SpProb>::remove(const Key k)
{
    auto ind = h(k);

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const SpProb*>(this)->mod(i);
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
            return true;
        }
    }
    return false;
}

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
