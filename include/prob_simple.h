#pragma once

#include "prob_base.h"

/* Fast LinProb using powers of 2 and bitmasking ******************************/
template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class FastLinProb : public ProbTraits<FastLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = FastLinProb<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename ProbTraits<This_t>::Key;
    using Data   = typename ProbTraits<This_t>::Data;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    static size_t next_power_of_two(size_t i)
    {
        size_t t = 2048;
        while (t < i) t <<= 1;
        return t;
    }

    FastLinProb(size_t cap = 0      , double /*size_constraint*/ = 1., size_t /*steps*/ = 0)
        : Base_t(next_power_of_two(cap) << 1, 1.)
    {
        thresh  = capacity * .6;
        bitmask = capacity - 1;
    }

    FastLinProb(const FastLinProb&) = delete;
    FastLinProb& operator=(const FastLinProb&) = delete;

    FastLinProb(FastLinProb&& rhs)  = default;
    FastLinProb& operator=(FastLinProb&& ) = default;

    inline size_t index(size_t i) const { return i & bitmask; }
    inline size_t mod(size_t i)   const { return i & bitmask; }

private:
    using Base_t::capacity;
    using Base_t::thresh;
    using Base_t::table;

    explicit FastLinProb(size_t cap, size_t lthresh, This_t* /* make unique */)
        : Base_t(cap, 1.), bitmask(cap-1)
    {
        thresh = lthresh;
    }

    inline void grow()
    {
        auto nsize  = capacity << 1;
        auto ntable = This_t(nsize, nsize*0.6, this);

        for (size_t i = 0; i <= bitmask; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        (*this) = std::move(ntable);
    }

    size_t bitmask;
};


template<class K, class D, class HF, class Conf>
class ProbTraits<FastLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = FastLinProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};



/* Using Classic Linear Probing to Fill a Table Densely ***********************/
template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class SpaceLinProb : public ProbTraits<SpaceLinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key    = typename ProbTraits<This_t>::Key;
    using Data   = typename ProbTraits<This_t>::Data;

    SpaceLinProb(size_t cap = 0, double size_constraint = 1.1, size_t /*steps*/=0)
        : Base_t(cap, size_constraint)
    {
        factor = double(capacity-300)/double(1ull << 32);
    }

    SpaceLinProb(const SpaceLinProb&) = delete;
    SpaceLinProb& operator=(const SpaceLinProb&) = delete;

    SpaceLinProb(SpaceLinProb&& rhs)  = default;
    SpaceLinProb& operator=(SpaceLinProb&& ) = default;

    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index(size_t i) const { return (bitmask & i)*factor; }
    inline size_t mod(size_t i)   const { return i; }

    using Base_t::alpha;
    using Base_t::n;
    using Base_t::capacity;

private:
    inline void grow()
    {
        auto ntable = This_t(n, alpha);

        for (size_t i = 0; i < capacity; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        (*this) = std::move(ntable);
    }

    using Base_t::table;

    double factor;
    static constexpr size_t bitmask = (1ull<<32)-1;
};


template<class K, class D, class HF, class Conf>
class ProbTraits<SpaceLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};
