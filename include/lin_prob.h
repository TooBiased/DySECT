#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "config.h"

template <class T>
class CuckooTraits;

template<class SpLinProb>
class LinProbBase
{
public:
    using Key      = typename CuckooTraits<SpLinProb>::Key;
    using Data     = typename CuckooTraits<SpLinProb>::Data;
    using FRet     = std::pair<bool, Data>;

private:
    using This_t         = LinProbBase<SpLinProb>;
    using Specialized_t  = typename CuckooTraits<SpLinProb>::Specialized_t;
    using Bucket_t       = Bucket<Key,Data,1>;
    using Cell_t   = std::pair<Key,Data>;
    using HashFct_t      = typename CuckooTraits<SpLinProb>::HashFct_t;
    using HistCount_t    = typename CuckooTraits<SpLinProb>::Config_t::HistCount_t;

    friend Specialized_t;

public:
    LinProbBase(size_t cap, size_t dis_steps)
        : n(0), capacity(cap), steps(dis_steps), hcounter(dis_steps)
    {
        if (cap) table = std::make_unique<Cell_t[]>(cap);
    }

    ~LinProbBase() = default;

    LinProbBase(const LinProbBase&) = delete;
    LinProbBase& operator=(const LinProbBase&) = delete;

    LinProbBase(LinProbBase&& rhs)
        : n(rhs.n), capacity(rhs.capacity), steps(rhs.steps),
          hcounter(rhs.steps), table(std::move(rhs.table))
    { }

    LinProbBase& operator=(LinProbBase&& rhs)
    {
        n = rhs.n;
        capacity = rhs.capacity;
        table = std::move(rhs.table);
        return *this;
    }

    void init(size_t cap)
    {
        capacity = cap;
        table = std::make_unique<Cell_t[]>(cap);
    }

    bool insert(Key k, Data d);
    bool insert(std::pair<Key,Data> t);
    FRet find  (Key k) const;
    bool remove(Key k);

    /*** some functions for easier load visualization (not in final product) **/
    std::pair<size_t, Bucket_t*> getTable(size_t)
    { return std::make_pair(0, nullptr); }
    void clearHist()
    { for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0; }

    /*** members that should become private at some point *********************/
    size_t     n;
    size_t     capacity;
    size_t     steps;
    HashFct_t  hasher;
    HistCount_t  hcounter;
    std::unique_ptr<Cell_t[]> table;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

private:
    /*** static polymorph functions *******************************************/
    inline void inc_n() { ++n; }
    inline void dec_n() { --n; }
};



/* IMPLEMENTATION *************************************************************/
template<class SpLinProb>
inline bool LinProbBase<SpLinProb>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class SpLinProb>
inline bool LinProbBase<SpLinProb>::insert(std::pair<Key, Data> t)
{
    auto ind = static_cast<const Specialized_t*>(this)->index(hasher(t.first));

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const SpLinProb*>(this)->mod(i);
        auto temp = table[ti];
        if ( temp.first == 0 )
        {
            table[ti] = t;
            hcounter.add(i - ind);
            static_cast<SpLinProb*>(this)->inc_n();
            return true;
        }
    }
    return false;
}

template<class SpLinProb>
inline typename LinProbBase<SpLinProb>::FRet
LinProbBase<SpLinProb>::find(Key k) const
{
    auto ind = static_cast<const Specialized_t*>(this)->index(hasher(k));

    for (size_t i = ind; ; ++i)
    {
        size_t ti = static_cast<const SpLinProb*>(this)->mod(i);
        auto temp = table[ti];

        if ( temp.first == 0 )
        {
            break;
        }
        else if ( temp.first == k )
        {
            return std::make_pair(true, temp.second);
        }
    }
    return std::make_pair(false, 0);
}

template<class SpLinProb>
inline bool LinProbBase<SpLinProb>::remove(Key k)
{
    /*
    auto ind = hash(k);
    size_t i = ind;
    while (true)
    {
        auto temp = table[i];
        if ( temp.first == 0 )
        {
            return false;
        }
        else if ( temp.first == k )
        {
            auto i2 = i+1;
            while (true)
            {
                auto temp2 = table[static_cast<SpLinProb*>(this)->mod(i2)];
                if (temp2.first == 0)
                    break;
                else if (before(hash(temp2.first), i))
                {
                    table[i & bitmask] = temp2;
                    i                  = i2;
                }
                i2++;
            }
            table[i & bitmask] = Cell_t(0, 0);
        }

        ++i;
    }
    */
    return false;
}


template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class FastLinProb : public CuckooTraits<FastLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = FastLinProb<K,D,HF,Conf>;
    using Base_t = typename CuckooTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename CuckooTraits<This_t>::Key;
    using Data   = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    FastLinProb(size_t cap = 0      , double /*size_constraint*/ = 1.1,
                size_t dis_steps = 0, size_t /*seed*/ = 0)
        : Base_t(0, dis_steps)
    {
        size_t c = 2048;
        while (c < cap) c <<= 1;
        c <<= 1;

        init(c);
        bitmask = c - 1;
        grow_thresh = c * 0.6;
    }

    FastLinProb(const FastLinProb&) = delete;
    FastLinProb& operator=(const FastLinProb&) = delete;

    FastLinProb(FastLinProb&& rhs)  = default;
    FastLinProb& operator=(FastLinProb&& ) = default;

    inline size_t index(size_t i) const { return i & bitmask; }
    inline size_t mod(size_t i)   const { return i & bitmask; }

private:
    using Base_t::capacity;
    using Base_t::init;
    using Base_t::steps;
    using Base_t::table;
    using Base_t::n;

    FastLinProb(size_t cap, size_t steps)
        : Base_t(cap, steps), bitmask(cap-1), grow_thresh(cap)
    { }

    inline void grow()
    {
        auto nsize  = capacity << 1;
        auto ntable = This_t(nsize, steps);

        for (size_t i = 0; i <= bitmask; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        std::swap(capacity, ntable.capacity);
        std::swap(bitmask , ntable.bitmask);
        std::swap(table   , ntable.table);
        grow_thresh = (capacity) * 0.6;
    }

    inline void inc_n() { ++n; if (n > grow_thresh) grow(); }

    size_t bitmask;
    size_t grow_thresh;
};


template<class K, class D, class HF, class Conf>
class CuckooTraits<FastLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = FastLinProb<K,D,HF,Conf>;
    using Base_t        = LinProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};


template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class SpaceLinProb : public CuckooTraits<SpaceLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t = typename CuckooTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename CuckooTraits<This_t>::Key;
    using Data   = typename CuckooTraits<This_t>::Data;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    SpaceLinProb(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t /*seed*/ = 0)
        : Base_t(cap*size_constraint, dis_steps)
    { }

    SpaceLinProb(const SpaceLinProb&) = delete;
    SpaceLinProb& operator=(const SpaceLinProb&) = delete;

    SpaceLinProb(SpaceLinProb&& rhs)  = default;
    SpaceLinProb& operator=(SpaceLinProb&& ) = default;

    inline size_t index(size_t i) const { return i % capacity; }
    inline size_t mod(size_t i)   const { return i % capacity; }

private:
    using Base_t::capacity;
};


template<class K, class D, class HF, class Conf>
class CuckooTraits<SpaceLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t        = LinProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};
