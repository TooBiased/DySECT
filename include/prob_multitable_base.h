#pragma once

#include "prob_simple.h"
#include "prob_robin.h"

namespace dysect
{

template<class K, class D, class HF, class Conf>
class multitable_linear
{
private:
    using This_t     = multitable_linear<K,D,HF,Conf>;
    using Subtable_t = prob_linear<K,D,HF,Conf>;
    using HashFct_t  = HF;

public:
    using key_type    = typename Subtable_t::key_type;
    using mapped_type = typename Subtable_t::mapped_type;

    using iterator    = typename Subtable_t::iterator;
    using const_iterator = typename Subtable_t::const_iterator;

    static constexpr size_t tl   = 256;
    static constexpr size_t bits = tl-1;

private:
    HashFct_t hasher;
    Subtable_t tables[tl];

public:
    multitable_linear(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            tables[i] = Subtable_t(cap/tl, size_constraint);
        }
    }

    inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
    {
        return tables[getInd(t.first)].insert(t);
    }

    inline iterator find(key_type k)
    {
        return tables[getInd(k)].find(k);
    }

    inline const_iterator find(key_type k) const
    {
        return tables[getInd(k)].find(k);
    }

    inline size_t erase(key_type k)
    {
        return tables[getInd(k)].erase(k);
    }

    inline int displacement(key_type k) const
    {
        return tables[getInd(k)].displacement(k);
    }


    inline iterator begin()              { return tables[0].begin(); }
    inline iterator end()                { return tables[0].end(); }
    inline const_iterator begin()  const { return tables[0].begin(); }
    inline const_iterator end()    const { return tables[0].end(); }
    inline const_iterator cbegin() const { return tables[0].cbegin(); }
    inline const_iterator cend()   const { return tables[0].cend(); }

private:
    inline size_t getInd(key_type k) const
    {
        return hasher(k) & bits;
    }

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(9); out << "f_cap" << " " << std::flush;
    }

    inline void print_init_data(std::ostream& out)
    {
        size_t cap = 0;
        for (size_t i = 0; i < tl; ++i)
            cap += tables[i].get_capacity();
        out.width(9); out << cap << " " << std::flush;
    }
};

template<class K, class D, class HF, class Conf>
class multitable_robin
{
private:
    using This_t     = multitable_robin<K,D,HF,Conf>;
    using Subtable_t = prob_robin<K,D,HF,Conf>;
    using HashFct_t  = HF;

public:
    using key_type    = typename Subtable_t::key_type;
    using mapped_type = typename Subtable_t::mapped_type;

    using iterator    = typename Subtable_t::iterator;
    using const_iterator = typename Subtable_t::const_iterator;

    static constexpr size_t tl   = 256;
    static constexpr size_t bits = tl-1;

private:
    HashFct_t hasher;
    Subtable_t tables[tl];

public:
    multitable_robin(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            tables[i] = Subtable_t(cap/tl, size_constraint);
        }
    }

    inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
    {
        return tables[getInd(t.first)].insert(t);
    }

    inline iterator find(key_type k)
    {
        return tables[getInd(k)].find(k);
    }

    inline size_t erase(key_type k)
    {
        return tables[getInd(k)].erase(k);
    }



    inline iterator begin()              { return tables[0].begin(); }
    inline iterator end()                { return tables[0].end(); }
    inline const_iterator begin()  const { return tables[0].begin(); }
    inline const_iterator end()    const { return tables[0].end(); }
    inline const_iterator cbegin() const { return tables[0].cbegin(); }
    inline const_iterator cend()   const { return tables[0].cend(); }

private:
    inline size_t getInd(key_type k)
    {
        return hasher(k) & bits;
    }

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(9); out << "f_cap" << " " << std::flush;
    }

    inline void print_init_data(std::ostream& out)
    {
        size_t cap = 0;
        for (size_t i = 0; i < tl; ++i)
            cap += tables[i].get_capacity();
        out.width(9); out << cap << " " << std::flush;
    }
};

};
