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
    using size_type   = typename Base_t::size_type;
    using key_type    = typename ProbTraits<This_t>::key_type;
    using mapped_type = typename ProbTraits<This_t>::mapped_type;


    FastLinProb(size_type cap = 0, double = 1., size_type = 0)
        : Base_t(next_power_of_two(cap) << 1, 1.)
    {
        thresh  = capacity * .6;
        bitmask = capacity - 1;
    }
    FastLinProb(const FastLinProb&) = delete;
    FastLinProb& operator=(const FastLinProb&) = delete;
    FastLinProb(FastLinProb&& rhs)  = default;
    FastLinProb& operator=(FastLinProb&& ) = default;

private:
    using Base_t::capacity;
    using Base_t::thresh;
    using Base_t::table;

    size_type bitmask;

    // Access Functions ********************************************************
    inline size_type index(size_type i) const { return i & bitmask; }
    inline size_type mod(size_type i)   const { return i & bitmask; }

    // Helper Function *********************************************************
    static size_type next_power_of_two(size_type i)
    {
        size_type t = 2048;
        while (t < i) t <<= 1;
        return t;
    }

    // Growing *****************************************************************
    explicit FastLinProb(size_type cap, size_type lthresh, This_t* )
        : Base_t(cap, 1.), bitmask(cap-1)
    {
        thresh = lthresh;
    }

    inline void grow()
    {
        auto nsize  = capacity << 1;
        auto ntable = This_t(nsize, nsize*0.6, this);

        for (size_type i = 0; i <= bitmask; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        (*this) = std::move(ntable);
    }
};




template<class K, class D, class HF, class Conf>
class ProbTraits<FastLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = FastLinProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using HashFct_t     = HF;
    using Config_t      = Conf;

    using key_type      = K;
    using mapped_type   = D;
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
    using size_type   = typename Base_t::size_type;
    using key_type    = typename ProbTraits<This_t>::key_type;
    using mapped_type = typename ProbTraits<This_t>::mapped_type;

    SpaceLinProb(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
        : Base_t(cap, size_constraint)
    {
        factor = double(capacity-300)/double(1ull << 32);
    }
    SpaceLinProb(const SpaceLinProb&) = delete;
    SpaceLinProb& operator=(const SpaceLinProb&) = delete;
    SpaceLinProb(SpaceLinProb&& rhs)  = default;
    SpaceLinProb& operator=(SpaceLinProb&& ) = default;

private:
    using Base_t::alpha;
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::table;

    double factor;

    static constexpr size_type bitmask = (1ull<<32)-1;

    // Access Functions ********************************************************
    inline size_type index(size_type i) const { return (bitmask & i)*factor; }
    inline size_type mod(size_type i)   const { return i; }

    // Growing *****************************************************************
    inline void grow()
    {
        auto ntable = This_t(n, alpha);

        for (size_type i = 0; i < capacity; ++i)
        {
            auto temp = table[i];
            if (temp.first)
            {
                ntable.insert(temp);
            }
        }

        (*this) = std::move(ntable);
    }
};


template<class K, class D, class HF, class Conf>
class ProbTraits<SpaceLinProb<K,D,HF,Conf> >
{
public:
    using Specialized_t = SpaceLinProb<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using HashFct_t     = HF;
    using Config_t      = Conf;

    using key_type      = K;
    using mapped_type   = D;
};










// *****************************************************************************
// Same as Above, but Growing Using in Place Migration *************************
// *****************************************************************************

template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class SpaceLinProbInPlace : public ProbTraits<SpaceLinProbInPlace<K,D,HF,Conf> >::Base_t
{
private:
    using This_t = SpaceLinProbInPlace<K,D,HF,Conf>;
    using Base_t = typename ProbTraits<This_t>::Base_t;

    friend Base_t;

public:
    using size_type   = typename Base_t::size_type;
    using key_type    = typename ProbTraits<This_t>::key_type;
    using mapped_type = typename ProbTraits<This_t>::mapped_type;
private:
    using value_intern = std::pair<key_type, mapped_type>;

    static constexpr size_type max_size = 16ull << 30;

public:
    SpaceLinProbInPlace(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
        : Base_t(0, size_constraint), bla(0)
    {
        // factor = double(capacity-300)/double(1ull << 32);

        value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
        if (temp) table = std::unique_ptr<value_intern[]>(temp);

        capacity = (cap) ? cap*alpha : 2048*alpha;
        thresh   = (cap) ? cap*beta  : 2048*beta;
        factor = double(capacity-300)/double(1ull << 32);

        std::fill(table.get(), table.get()+capacity, value_intern());
    }
    SpaceLinProbInPlace(const SpaceLinProbInPlace&) = delete;
    SpaceLinProbInPlace& operator=(const SpaceLinProbInPlace&) = delete;
    SpaceLinProbInPlace(SpaceLinProbInPlace&& rhs)  = default;
    SpaceLinProbInPlace& operator=(SpaceLinProbInPlace&& ) = default;

private:
    using Base_t::alpha;
    using Base_t::beta;
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::thresh;
    using Base_t::table;

    double    factor;
    size_type bla;

    static constexpr size_type bitmask = (1ull<<32)-1;

    using Base_t::h;

    // Access Functions ********************************************************
    inline size_type index(size_type i) const { return (bitmask & i)*factor; }
    inline size_type mod(size_type i)   const { return i; }

public:
    using Base_t::insert;

private:
    // Growing *****************************************************************
    inline void grow()
    {
        // auto ntable = This_t(n, alpha);

        size_type ncap    = n*alpha;
        size_type nthresh = n*beta;
        double    nfactor = double(ncap-300)/double(1ull << 32);

        std::fill(table.get()+capacity, table.get()+ncap, value_intern());

        capacity = ncap;
        thresh   = nthresh;
        factor   = nfactor;
        n = 0;

        std::vector<value_intern> buffer;

        for (int i = capacity - 1; i >= 0; --i)
        {
            auto temp = table[i];

            if (temp.first)
            {
                table[i] = value_intern();
                auto ind = h(temp.first);
                if (ind >= size_t(i))
                    insert(temp);
                else
                    buffer.push_back(temp);
            }
            else if (! buffer.empty())
            {
                bla = std::max(bla, buffer.size());
                for (auto it = buffer.begin(); it != buffer.end(); it++)
                {
                    insert(*it);
                }
                buffer.clear();
            }
        }

        if (! buffer.empty())
        {
            bla = std::max(bla, buffer.size());
            for (auto it = buffer.begin(); it != buffer.end(); it++)
            {
                insert(*it);
            }
            buffer.clear();
        }
    }

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(6); out  << "busize" << " ";
        Base_t::print_init_header(out);
    }

    inline        void print_init_data  (std::ostream& out)
    {
        out.width(6);  out << bla << " ";
        Base_t::print_init_data(out);
    }
};


template<class K, class D, class HF, class Conf>
class ProbTraits<SpaceLinProbInPlace<K,D,HF,Conf> >
{
public:
    using Specialized_t = SpaceLinProbInPlace<K,D,HF,Conf>;
    using Base_t        = ProbBase<Specialized_t>;
    using HashFct_t     = HF;
    using Config_t      = Conf;

    using key_type      = K;
    using mapped_type   = D;
};
