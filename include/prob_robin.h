#pragma once

#include "prob_base.h"

template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class RobinProb : public ProbTraits<RobinProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = RobinProb<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using Key      = typename ProbTraits<This_t>::Key;
    using Data     = typename ProbTraits<This_t>::Data;
    using iterator = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

    RobinProb(size_t cap = 0      , double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : Base_t(cap, size_constraint),
          pdistance(0)
    {
        factor = double(capacity-300)/double(1ull << 32);
    }

    RobinProb(const RobinProb&) = delete;
    RobinProb& operator=(const RobinProb&) = delete;

    RobinProb(RobinProb&& rhs)  = default;
    RobinProb& operator=(RobinProb&& ) = default;

private:
    using Base_t::alpha;
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::hasher;
    using Base_t::table;
    using Base_t::h;
    using Base_t::inc_n;
    using Base_t::make_iterator;
    using Base_t::make_citerator;

    static constexpr size_t bitmask = (1ull << 32) - 1;
    size_t pdistance;
    double factor;

public:
    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index (size_t i) const
    { return double(bitmask & i) * factor; }
    inline double dindex(size_t i) const
    { return double(bitmask & i) * factor; }
    inline size_t mod(size_t i)    const
    { return i; }

    //specialized functions because of Robin Hood Hashing
    inline std::pair<iterator, bool> insert(const Key& k, const Data& d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator, bool> insert(const std::pair<Key, Data>& t)
    {
        // using doubles makes the element order independent from the capacity
        // thus growing gets even easier
        double ind     = dindex(hasher(t.first));
        auto   current = t;

        for (size_t i = ind; ; ++i)
        {
            auto ti    = mod(i);
            auto temp  = table[ti];

            if ( temp.first == t.first )
            {
                return std::make_pair(make_iterator(&table[ti]), false);
            }
            if ( temp.first == 0 )
            {
                if (i == capacity - 1)
                    return std::make_pair(Base_t::end(), false);
                table[ti] = current;
                inc_n();
                pdistance = std::max<size_t>(pdistance, i-size_t(ind));
                return std::make_pair(make_iterator(&table[ti]), true);
            }
            double tind = dindex(hasher(temp.first));
            if ( tind > ind )
            {
                std::swap(table[ti], current);
                pdistance = std::max<int>(pdistance, i-size_t(ind));
                ind = tind;
            }
        }
        return std::make_pair(Base_t::end(), false);
    }

    inline iterator find(const Key& k)
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
                return make_iterator(&table[ti]);
            }
        }
        return Base_t::end();
    }

    inline const_iterator find(const Key& k) const
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
                return make_citerator(&table[ti]);
            }
        }
        return Base_t::cend();
    }

    inline size_t erase(const Key& k)
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
                return 1;
            }
        }
        return 0;
    }

private:
    inline void grow()
    {
        auto ntable = This_t(n, alpha);

        size_t tn = n;
        size_t tdistance = ntable.migrate(*this);

        (*this) = std::move(ntable);

        n         = tn;
        pdistance = tdistance;
    }

    inline size_t migrate(This_t& source)
    {
        size_t distance   = 0;
        size_t target_pos = 0;
        for (size_t i = 0; i < source.capacity; ++i)
        {
            auto current = source.table[i];
            if (!current.first) continue;
            auto hash = h(current.first);
            if (target_pos > hash)
                distance   = std::max<size_t>(distance, target_pos - hash);
            else
                target_pos = hash;
            table[target_pos++] = current;
        }
        return distance;
    }

    inline void propagate_remove(const size_t hole)
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

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(5); out  << "pdis" << " ";
        Base_t::print_init_header(out);
    }

    inline        void print_init_data  (std::ostream& out)
    {
        out.width(5);  out << pdistance << " ";
        Base_t::print_init_data(out);
    }
};


template<class K, class D, class HF, class Conf>
class ProbTraits<RobinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = RobinProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};
