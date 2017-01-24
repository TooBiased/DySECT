#pragma once

#include "prob_base.h"


template <size_t ns>
class uiType;

template <>
class uiType<64>
{ using type = size_t; };

template <>
class uiType<32>
{ using type = uint; };

template <>
class uiType<16>
{ using type = ushort; };

template <size_t ns = 64>
class AugmentMap
{
public:
    using size_type  = typename uiType<ns>::type;
    size_type bitmap;

    AugmentMap() : bitmap(0) {}

    static constexpr size_t nh_size = ns;

    inline bool check(size_t i)
    {
        return (1ull << i) & bitmap;
    }

    inline void set(size_t i)
    {
        bitmap |= (1ull << i);
    }

    inline void unset(size_t i)
    {
        bitmap &= ~(1ull << i);
    }
};


template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class HopsProb : public ProbTraits<HopsProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = HopsProb<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;
    using Augment_t = AugmentMap<64>;

    friend Base_t;

public:
    using Key    = typename ProbTraits<This_t>::Key;
    using Data   = typename ProbTraits<This_t>::Data;

    HopsProb(size_t cap = 0      , double size_constraint = 1.1,
            size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : Base_t(cap, size_constraint)
    {
        factor = double(capacity-Augment_t::nh_size)/double(1ull << 32);
        nh_data = std::make_unique<Augment_t[]>(capacity-Augment_t::nh_size+1);
    }

    HopsProb(const HopsProb&) = delete;
    HopsProb& operator=(const HopsProb&) = delete;

    HopsProb(HopsProb&& rhs)  = default;
    HopsProb& operator=(HopsProb&& ) = default;

    /* SHOULD CHANGE THIS TO THE MULTIPLY BY DOUBLE FACTOR VARIANT */
    inline size_t index (size_t i) const
    { return double(bitmask & i) * factor; }
    inline size_t mod(size_t i)    const
    { return i; } //i%capacity; }

private:
    inline void grow()
    {
        auto ntable = This_t(n, alpha);

        for (size_t i = 0; i < capacity; ++i)
        {
            auto current  = table[i];
            ntable.insert(current);
        }

        (*this) = std::move(ntable);
    }

    using Base_t::alpha;
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::table;

    using Base_t::h;
    using Base_t::inc_n;

    static constexpr size_t bitmask = (1ull << 32) - 1;
    std::unique_ptr<Augment_t[]> nh_data;
    double factor;

public:
    //specialized functions because of Hop Hood Hashing
    inline bool insert(Key k, Data d)
    {
        return insert(std::make_pair(k,d));
    }

    inline bool insert(std::pair<Key, Data> t)
    {
        // we first have to check if t.first is already present
        size_t ind     = h(t.first);
        size_t bits  = nh_data[ind].bitmap;

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits & 1)) continue;
            else
            {
                auto temp = table[i];
                if ( temp.first == t.first )
                {
                    return false;
                }
            }
        }

        for (size_t i = ind; ; ++i)
        {
            auto temp  = table[i];
            if ( temp.first == 0)
            {
                size_t ti = i;
                if (ti >= ind + Augment_t::nh_size)
                {
                    bool successful;
                    std::tie(successful, ti) = move_gap(i, ind+Augment_t::nh_size);
                    if (!successful) break;
                }
                table[ti] = t;
                inc_n();
                nh_data[ind].set(ti - ind);
                return true;
            }
        }
        return false;
    }

    inline std::pair<bool, size_t> move_gap(size_t pos, size_t goal)
    {
        for (size_t i = pos - Augment_t::nh_size + 1; i < pos; ++i)
        {
            auto current = table[i];
            auto ind = h(current.first);
            if (ind + Augment_t::nh_size > pos)
            {
                nh_data[ind].unset(i - ind);
                nh_data[ind].set(pos - ind);
                table[pos] = current;
                table[i]   = std::make_pair(0,0);
                if (i < goal) return std::make_pair(true, i);
                else return move_gap(i, goal);
            }
        }
        return std::make_pair(false, pos);
    }

    inline typename Base_t::FRet find(Key k) const
    {
        auto ind = h(k);
        size_t bits  = nh_data[ind].bitmap;

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits & 1)) continue;
            else
            {
                auto temp = table[i];
                if ( temp.first == k )
                {
                    return std::make_pair(true, temp.second);
                }
            }
        }
        return std::make_pair(false, 0);
    }

    inline bool remove(Key k)
    {
        auto ind = h(k);
        auto nh_ind = nh_data[ind];

        for (size_t i = 0; i < Augment_t::nh_size; ++i)
        {
            if (!nh_ind.check(i)) continue;
            auto temp = table[ind+i];
            if ( temp.first == k )
            {
                nh_data[ind].unset(i);
                table[ind+i] = std::make_pair(0,0);
                return true;
            }
        }
        return false;
    }

};


template<class K, class D, class HF, class Conf>
class ProbTraits<HopsProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = HopsProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using Key           = K;
    using Data          = D;
    using HashFct_t     = HF;
    using Config_t      = Conf;
};
