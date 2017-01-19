#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "config.h"
#include "bucket.h"

template <class T>
class LinProbTraits;

template<class SpLinProb>
class LinProbBase
{
public:
    using Key      = typename LinProbTraits<SpLinProb>::Key;
    using Data     = typename LinProbTraits<SpLinProb>::Data;
    using FRet     = std::pair<bool, Data>;

private:
    using This_t         = LinProbBase<SpLinProb>;
    using Specialized_t  = typename LinProbTraits<SpLinProb>::Specialized_t;
    using Bucket_t       = Bucket<Key,Data,1>;
    using Cell_t   = std::pair<Key,Data>;
    using HashFct_t      = typename LinProbTraits<SpLinProb>::HashFct_t;
    using HistCount_t    = typename LinProbTraits<SpLinProb>::Config_t::HistCount_t;

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
    bool remove(Key k);// { return false; }

    /*** some functions for easier load visualization (not in final product) **/
    std::pair<size_t, Bucket_t*> getTable(size_t)
    { return std::make_pair(0, nullptr); }
    void clearHist()
    { for (size_t i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0; }
    inline static void print_init_header(std::ostream& out)
        { out.width(9); out << "f_cap"  << std::flush;}
    inline void print_init_data  (std::ostream& out)
        { out.width(9); out << capacity << std::flush;}

    /*** members that should become private at some point *********************/
    size_t     n;
    size_t     capacity;
    size_t     steps;
    HashFct_t  hasher;
    HistCount_t  hcounter;
    std::unique_ptr<Cell_t[]> table;

private:
    /*** static polymorph functions *******************************************/
    inline size_t h(Key k) const
    { return static_cast<const Specialized_t*>(this)->index(hasher(k)); }
    void propagate_remove(size_t origin);
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
    auto ind = h(t.first);

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
    auto ind = h(k);

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
    auto ind = h(k);

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
            dec_n();
            table[ti] = std::make_pair(0,0);
            propagate_remove(ti);
            return true;
        }
    }
    return false;
}

template<class SpLinProb>
inline void LinProbBase<SpLinProb>::propagate_remove(size_t origin)
{
    size_t tempn = n;
    n = 0;
    for (size_t i = origin+1; ; ++i)
    {
        size_t ti = static_cast<const SpLinProb*>(this)->mod(i);
        auto temp = table[ti];

        if (temp.first == 0) break;

        table[ti] = std::make_pair(0,0);
        insert(temp);
    }
    n = tempn;
}

template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class FastLinProb : public LinProbTraits<FastLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = FastLinProb<K,D,HF,Conf>;
    using Base_t = typename LinProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename LinProbTraits<This_t>::Key;
    using Data   = typename LinProbTraits<This_t>::Data;

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
class LinProbTraits<FastLinProb<K,D,HF,Conf> >
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
class SpaceLinProb : public LinProbTraits<SpaceLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t = typename LinProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename LinProbTraits<This_t>::Key;
    using Data   = typename LinProbTraits<This_t>::Data;

    SpaceLinProb(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t /*seed*/ = 0)
        : Base_t(0, dis_steps), alpha(size_constraint), beta((alpha+1.)/2.)
    {
        size_t c = std::max(cap, size_t(2048));
        init(c*alpha);
        thresh = c*beta;
    }

    SpaceLinProb(const SpaceLinProb&) = delete;
    SpaceLinProb& operator=(const SpaceLinProb&) = delete;

    SpaceLinProb(SpaceLinProb&& rhs)  = default;
    SpaceLinProb& operator=(SpaceLinProb&& ) = default;

    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index(size_t i) const { return i % capacity; }
    inline size_t mod(size_t i)   const { return i % capacity; }

private:
    inline void grow()
    {
        //auto nsize  = n*alpha;
        auto ntable = This_t(n, alpha, steps);

        for (size_t i = 0; i < capacity; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        std::swap(capacity   , ntable.capacity);
        std::swap(table      , ntable.table);
        std::swap(thresh     , ntable.thresh);
    }

    inline void inc_n() { ++n; if (n > thresh) grow(); }

    using Base_t::capacity;
    using Base_t::n;
    using Base_t::table;
    using Base_t::steps;
    double alpha;
    double beta;
    size_t thresh;

    using Base_t::init;
};


template<class K, class D, class HF, class Conf>
class LinProbTraits<SpaceLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t        = LinProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};


template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class RobinProb : public LinProbTraits<RobinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = RobinProb<K,D,HF,Conf>;
    using Base_t = typename LinProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename LinProbTraits<This_t>::Key;
    using Data   = typename LinProbTraits<This_t>::Data;

    RobinProb(size_t cap = 0      , double size_constraint = 1.1,
                 size_t dis_steps = 0, size_t /*seed*/ = 0)
        : Base_t(0, dis_steps), alpha(size_constraint), beta((alpha+1.)/2.),
          pdistance(0)
    {
        size_t c = std::max(cap, size_t(2048));
        thresh = c*beta;

        c = size_t(double(c)*alpha);
        init(c);
        factor = double(c)/double(1ull << 32);
    }

    RobinProb(const RobinProb&) = delete;
    RobinProb& operator=(const RobinProb&) = delete;

    RobinProb(RobinProb&& rhs)  = default;
    RobinProb& operator=(RobinProb&& ) = default;

    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index (size_t i) const { return i%capacity; }
    inline double dindex(size_t i) const { return i%capacity; }
    inline size_t mod(size_t i)   const  { return i%capacity; }

private:
    inline void grow()
    {
        auto ntable = This_t(n, alpha, steps);

        for (size_t i = 0; i < capacity; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.Base_t::insert(temp);
            }
        }

        std::swap(capacity, ntable.capacity);
        std::swap(table   , ntable.table);
        std::swap(thresh  , ntable.thresh);
    }

    inline void inc_n() { ++n; if (n > thresh) grow(); }

    using Base_t::capacity;
    using Base_t::n;
    using Base_t::table;
    using Base_t::steps;
    using Base_t::hasher;
    using Base_t::h;
    double alpha;
    double beta;
    double factor;
    size_t thresh;
    size_t pdistance;

    using Base_t::init;
public:
    //specialized functions because of Robin Hood Hashing
    inline bool insert(Key k, Data d)
    {
        return insert(std::make_pair(k,d));
    }

    inline bool insert(std::pair<Key, Data> t)
    {
        size_t ind     = h(t.first);
        auto   current = t;

        for (size_t i = ind; ; ++i)
        {
            size_t ti  = mod(i);
            auto temp  = table[ti];
            if ( temp.first == 0 )
            {
                table[ti] = current;
                inc_n();
                pdistance = std::max(pdistance, size_t(i-ind));
                return true;
            }
            auto tind = index(hasher(temp.first));
            if ( tind > ind )
            {
                std::swap(table[ti], current);
                pdistance = std::max(pdistance, size_t(i-ind));
                ind = tind;
            }
        }
        return false;
    }

    inline typename Base_t::FRet find(Key k) const
    {
        auto ind = h(k);

        for (size_t i = ind; i < ind+pdistance; ++i)
        {
            size_t ti = mod(i);
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

    inline bool remove(Key k)
    {
        auto ind = h(k);

        for (size_t i = ind; i < ind+pdistance ; ++i)
        {
            size_t ti = mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 )
            {
                break;
            }
            else if ( temp.first == k )
            {
                Base_t::dec_n();
                propagate_remove(ti);
                return true;
            }
        }
        return false;
    }

    inline void propagate_remove(size_t hole)
    {
        size_t thole = hole;
        for (size_t i = hole+1; ; ++i)
        {
            auto ti   = mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 ) break;
            auto tind = h(temp.first);
            if ( tind >= ti)       break;

            table[thole] = temp;
            thole = ti;
        }
        table[thole] = std::make_pair(0,0);
    }

    inline void print_init_data  (std::ostream& out)
    {
        out.width(9); out << capacity;
        out.width(15); out << pdistance << std::flush;
    }
};


template<class K, class D, class HF, class Conf>
class LinProbTraits<RobinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = RobinProb<K,D,HF,Conf>;
    using Base_t        = LinProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};
