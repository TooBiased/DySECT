#pragma once

#include "prob_base.h"

template <class AugmentData>
class AugmentDataAccessor
{
public:
    static constexpr size_t bitmask = AugmentData::bitmask;
    static constexpr size_t nh_size = AugmentData::nh_size;

    AugmentDataAccessor(size_t* data) : data(*data) { }

    size_t& data;

    inline size_t getNHood() const
    {
        return data & bitmask;
    }
    inline void set(size_t i)
    {
        data |= (1ull << i);
    }
    inline void unset(size_t i)
    {
        data &= ~(1ull << i);
    }
};

template <size_t ns>
class AugmentData
{
    using uchar = uint8_t;
public:
    static_assert(ns >   0, "Augment Data cannot handle neighborhood size 0!");
    static_assert(ns <= 64, "Augment Data cannot handle neighborhoods >64!");

    static constexpr size_t bitmask = (ns == 64) ? ~size_t(0) : (1ull << ns) -1;
    static constexpr size_t n_bytes = (ns>>3)+((ns&7) ? 1:0);
    static constexpr size_t nh_size = ns;
    using This_t     = AugmentData<ns>;
    using Accessor_t = AugmentDataAccessor<This_t>;

    AugmentData(size_t capacity)
    {
        auto glob_n_bytes = capacity*n_bytes+8-n_bytes;
        //data = (uchar*) malloc(glob_n_bytes);
        data = std::make_unique<uchar[]>(glob_n_bytes);
        //std::fill(data, data+glob_n_bytes, 0);
    }

    AugmentData(const AugmentData&) = delete;
    AugmentData& operator=(const AugmentData&) = delete;

    AugmentData(AugmentData&& rhs) : data(nullptr)
    {   std::swap(data, rhs.data); }
    AugmentData& operator=(AugmentData&& rhs)
    {   std::swap(data, rhs.data); return *this; }

    inline Accessor_t getAcc(size_t index)
    { return Accessor_t(reinterpret_cast<size_t*>(data.get() + index*n_bytes)); }

    inline const Accessor_t getCAcc(size_t index) const
    { return Accessor_t(reinterpret_cast<size_t*>(data.get() + index*n_bytes)); }

    inline size_t getNHood(size_t index) const
    {
        size_t temp = *reinterpret_cast<const size_t*>(data.get() + index*n_bytes);
        return temp & bitmask;
    }

    std::unique_ptr<uchar[]> data;
};


template <class K, class D, class HF = std::hash<K>,
          class Conf = HopscotchConfig<> >
class HopsProb : public ProbTraits<HopsProb<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = HopsProb<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;

    static constexpr size_t nh_size = Conf::NeighborSize;
    using AugData_t = AugmentData<nh_size>;

    friend Base_t;

public:
    using key_type      = typename ProbTraits<This_t>::key_type;
    using mapped_type     = typename ProbTraits<This_t>::mapped_type;
    using iterator = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

    HopsProb(size_t cap = 0      , double size_constraint = 1.1,
            size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : Base_t(cap, size_constraint), nh_data(capacity-nh_size+1)
    {
        factor = double(capacity-nh_size)/double(1ull << 32);
    }

    HopsProb(const HopsProb&) = delete;
    HopsProb& operator=(const HopsProb&) = delete;

    HopsProb(HopsProb&& rhs)  = default;
    HopsProb& operator=(HopsProb&& ) = default;

private:
    using Base_t::alpha;
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::table;

    double    factor;
    AugData_t nh_data;

    static constexpr size_t bitmask = (1ull << 32) - 1;

    using Base_t::h;
    using Base_t::inc_n;
    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    //specialized functions because of Hops Hashing
    inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>& t)
    {
        // we first have to check if t.first is already present
        size_t ind  = h(t.first);
        auto   aug  = nh_data.getAcc(ind);
        size_t bits = aug.getNHood();

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits & 1)) continue;
            else
            {
                auto temp = table[i];
                if ( temp.first == t.first )
                {
                    return std::make_pair(make_iterator(&table[i]), false);
                }
            }
        }

        for (size_t i = ind; ; ++i)
        {
            auto temp  = table[i];
            if ( temp.first == 0)
            {
                size_t ti = i;
                if (ti >= ind + nh_size)
                {
                    bool successful;
                    std::tie(successful, ti) = move_gap(i, ind+nh_size);
                    if (!successful) break;
                }
                table[ti] = t;
                aug.set(ti-ind);
                inc_n();
                return std::make_pair(make_iterator(&table[ti]), true);
            }
        }
        return std::make_pair(Base_t::end(), false);
    }

    inline iterator find(const key_type& k)
    {
        auto ind = h(k);
        size_t bits = nh_data.getNHood(ind);

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits & 1)) continue;
            auto temp = table[i];
            if ( temp.first == k )
            {
                return make_iterator(&table[i]);
            }
        }
        return Base_t::end();
    }

    inline const_iterator find(const key_type& k) const
    {
        auto ind = h(k);
        size_t bits = nh_data.getNHood(ind);

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits & 1)) continue;
            auto temp = table[i];
            if ( temp.first == k )
            {
                return make_citerator(&table[i]);
            }
        }
        return Base_t::cend();
    }

    inline size_t erase(const key_type& k)
    {
        auto ind = h(k);
        size_t bits = nh_data.getNHood(ind);

        for (size_t i = ind; bits; ++i, bits>>=1)
        {
            if (!(bits&1)) continue;
            auto tempk = table[i].first;
            if ( tempk == k )
            {
                nh_data.getAcc(ind).unset(i-ind);
                table[i] = std::make_pair(0,0);
                return 1;
            }
        }
        return 0;
    }

private:
    inline size_t index(size_t i) const
    { return double(bitmask & i) * factor; }
    inline size_t mod  (size_t i) const
    { return i; }

    inline void grow()
    {
        auto ntable = This_t(n, alpha);

        for (size_t i = 0; i < capacity; ++i)
        {
            auto current  = table[i];
            if (current.first)
                ntable.insert(current);
        }

        (*this) = std::move(ntable);
    }

    inline std::pair<bool, size_t> move_gap(const size_t pos, const size_t goal)
    {
        for (size_t i = pos - nh_size + 1; i < pos; ++i)
        {
            auto current = table[i];
            auto ind = h(current.first);
            if (ind + nh_size > pos)
            {
                auto aug = nh_data.getAcc(ind);
                aug.unset(i  -ind);
                aug.set  (pos-ind);

                table[pos] = current;
                table[i]   = std::make_pair(0,0);
                if (i < goal) return std::make_pair(true, i);
                else return move_gap(i, goal);
            }
        }
        return std::make_pair(false, pos);
    }

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(5); out << "nghb" << " ";
        Base_t::print_init_header(out);
    }

    inline void print_init_data(std::ostream& out)
    {
        out.width(5); out << nh_size << " ";
        Base_t::print_init_data(out);
    }
};


template<class K, class D, class HF, class Conf>
class ProbTraits<HopsProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = HopsProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using HashFct_t     = HF;
    using Config_t      = Conf;

    using key_type      = K;
    using mapped_type   = D;
};
