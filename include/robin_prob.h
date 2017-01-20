#pragma once

#include "lin_prob.h"

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
        thresh = (c-300)*beta;

        c = size_t(double(c)*alpha);
        init(c);
        factor = double(c-300)/double(1ull << 32);
    }

    RobinProb(const RobinProb&) = delete;
    RobinProb& operator=(const RobinProb&) = delete;

    RobinProb(RobinProb&& rhs)  = default;
    RobinProb& operator=(RobinProb&& ) = default;

    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index (size_t i) const
    {
        return double(bitmask & i) * factor;
    }
    inline double dindex(size_t i) const
    {
        return double(bitmask & i) * factor;
    }
    inline size_t mod(size_t i)   const  { return i; } //i%capacity; }

private:
    inline void grow()
    {
        auto ntable = This_t(n, alpha, steps);

        ntable.migrate(*this);

        std::swap(capacity , ntable.capacity);
        std::swap(table    , ntable.table);
        std::swap(thresh   , ntable.thresh);
        //std::swap(pdistance, ntable.pdistance);
        std::swap(factor   , ntable.factor);
    }

    inline void migrate(This_t& source)
    {
        size_t target_pos = 0;
        for (size_t i = 0; i < source.capacity; ++i)
        {
            auto current = source.table[i];
            if (!current.first) continue;
            auto hash = h(current.first);
            target_pos = std::max(hash, target_pos);
            table[target_pos++] = current;
        }

    }

    inline void inc_n() { ++n; if (n > thresh) grow(); }

    using Base_t::capacity;
    using Base_t::n;
    using Base_t::table;
    using Base_t::steps;
    using Base_t::hasher;
    using Base_t::h;
    static constexpr size_t bitmask = (1ull << 32) - 1;
    double alpha;
    double beta;
    size_t thresh;
    size_t pdistance;
    double factor;

    using Base_t::init;
public:
    //specialized functions because of Robin Hood Hashing
    inline bool insert(Key k, Data d)
    {
        return insert(std::make_pair(k,d));
    }

    inline bool insert(std::pair<Key, Data> t)
    {
        // using doubles makes the element order independent from the capacity
        // thus growing gets even easier
        double ind     = dindex(hasher(t.first));
        auto   current = t;

        for (size_t i = ind; ; ++i)
        {
            auto ti    = mod(i);
            auto temp  = table[ti];
            if ( temp.first == 0 )
            {
                if (i == capacity - 1) return false;
                table[ti] = current;
                inc_n();
                pdistance = std::max<int>(pdistance, i-ind);
                return true;
            }
            double tind = dindex(hasher(temp.first));
            if ( tind > ind )
            {
                std::swap(table[ti], current);
                pdistance = std::max<int>(pdistance, i-ind);
                ind = tind;
            }
        }
        return false;
    }

    inline typename Base_t::FRet find(Key k) const
    {
        auto ind = h(k);

        for (size_t i = ind; i <= ind+pdistance; ++i)
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

        for (size_t i = ind; i <= ind+pdistance ; ++i)
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

    inline static void print_init_header(std::ostream& out)
    {
        out.width(9); out  << "f_cap";
        out.width(5); out  << "pdis" << std::flush;
    }

    inline        void print_init_data  (std::ostream& out)
    {
        out.width(9);  out << capacity;
        out.width(5);  out << pdistance << std::flush;
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
