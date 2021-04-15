#pragma once

#include "utils/output.hpp"

#include "prob_simple.h"
#include "prob_robin.h"
#include "prob_quadratic.h"

namespace otm = utils_tm::out_tm;


// THESE TABLES USE THE LOWERMOST FEW BITS FOR THEIR TABLE ADDRESSING:
//  -> SUBTABLES SHOULD USE LINEAR MAPPING (MOSTLY USE UPPER BITS)



namespace dysect
{

template<class K, class D, class HF, class Conf>
class multitable_linear
{
private:
    using this_type           = multitable_linear<K,D,HF,Conf>;
    using subtable_type       = prob_linear<K,D,HF,Conf>;
    using hash_function_type  = HF;

public:
    using key_type       = typename subtable_type::key_type;
    using mapped_type    = typename subtable_type::mapped_type;

    using iterator       = typename subtable_type::iterator;
    using const_iterator = typename subtable_type::const_iterator;

    static constexpr size_t tl   = 256;
    static constexpr size_t bits = tl-1;

private:
    hash_function_type hasher;
    subtable_type tables[tl];

public:
    multitable_linear(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            tables[i] = subtable_type(cap/tl, size_constraint);
        }
    }

    inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
    {
        return tables[get_table(t.first)].insert(t);
    }

    inline iterator find(key_type k)
    {
        return tables[get_table(k)].find(k);
    }

    inline const_iterator find(key_type k) const
    {
        return tables[get_table(k)].find(k);
    }

    inline size_t erase(key_type k)
    {
        return tables[get_table(k)].erase(k);
    }

    inline int displacement(key_type k) const
    {
        return tables[get_table(k)].displacement(k);
    }


    inline iterator begin()              { return tables[0].begin(); }
    inline iterator end()                { return tables[0].end(); }
    inline const_iterator begin()  const { return tables[0].begin(); }
    inline const_iterator end()    const { return tables[0].end(); }
    inline const_iterator cbegin() const { return tables[0].cbegin(); }
    inline const_iterator cend()   const { return tables[0].cend(); }

private:
    inline size_t get_table(key_type k) const
    {
        return hasher(k) & bits;
    }

public:
    inline static void print_init_header(otm::output_type& out)
    {
        out << otm::width(10) << "f_cap";
    }

    inline void print_init_data(otm::output_type& out)
    {
        size_t cap = 0;
        for (size_t i = 0; i < tl; ++i)
            cap += tables[i].get_capacity();

        out << otm::width(10) << cap;
    }
};

template<class K, class D, class HF, class Conf>
class multitable_robin
{
private:
    using this_type          = multitable_robin<K,D,HF,Conf>;
    using subtable_type      = prob_robin<K,D,HF,Conf>;
    using hash_function_type = HF;

public:
    using key_type       = typename subtable_type::key_type;
    using mapped_type    = typename subtable_type::mapped_type;

    using iterator       = typename subtable_type::iterator;
    using const_iterator = typename subtable_type::const_iterator;

    static constexpr size_t tl   = 256;
    static constexpr size_t bits = tl-1;

private:
    hash_function_type hasher;
    subtable_type      tables[tl];

public:
    multitable_robin(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            tables[i] = subtable_type(cap/tl, size_constraint);
        }
    }

    inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
    {
        return tables[get_table(t.first)].insert(t);
    }

    inline iterator find(key_type k)
    {
        return tables[get_table(k)].find(k);
    }

    inline iterator find(key_type k) const
    {
        return tables[get_table(k)].find(k);
    }

    inline size_t erase(key_type k)
    {
        return tables[get_table(k)].erase(k);
    }

    inline int displacement(key_type k) const
    {
        return tables[get_table(k)].displacement(k);
    }

    inline iterator begin()              { return tables[0].begin(); }
    inline iterator end()                { return tables[0].end(); }
    inline const_iterator begin()  const { return tables[0].begin(); }
    inline const_iterator end()    const { return tables[0].end(); }
    inline const_iterator cbegin() const { return tables[0].cbegin(); }
    inline const_iterator cend()   const { return tables[0].cend(); }

private:
    inline size_t get_table(key_type k) const
    {
        return hasher(k) & bits;
    }

public:
    inline static void print_init_header(otm::output_type& out)
    {
        out << otm::width(10) << "f_cap";
    }

    inline void print_init_data(otm::output_type& out)
    {
        size_t cap = 0;
        for (size_t i = 0; i < tl; ++i)
            cap += tables[i].get_capacity();

        out << otm::width(10) << cap;
    }
};




template<class K, class D, class HF, class Conf>
class multitable_quadratic
{
private:
    using this_type          = multitable_quadratic<K,D,HF,Conf>;
    using subtable_type      = prob_quadratic<K,D,HF,Conf>;
    using hash_function_type = HF;

public:
    using key_type       = typename subtable_type::key_type;
    using mapped_type    = typename subtable_type::mapped_type;

    using iterator       = typename subtable_type::iterator;
    using const_iterator = typename subtable_type::const_iterator;

    static constexpr size_t tl   = 256;
    static constexpr size_t bits = tl-1;

private:
    hash_function_type hasher;
    subtable_type      tables[tl];

public:
    multitable_quadratic(size_t cap   = 0, double size_constraint = 1.1, size_t /**/ = 0, size_t /**/ = 0)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            tables[i] = subtable_type(cap/tl, size_constraint);
        }
    }

    inline std::pair<iterator,bool> insert(key_type k, mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline std::pair<iterator,bool> insert(std::pair<key_type,mapped_type> t)
    {
        return tables[get_table(t.first)].insert(t);
    }

    inline iterator find(key_type k)
    {
        return tables[get_table(k)].find(k);
    }

    inline iterator find(key_type k) const
    {
        return tables[get_table(k)].find(k);
    }

    inline size_t erase(key_type k)
    {
        return tables[get_table(k)].erase(k);
    }

    inline int displacement(key_type k) const
    {
        return tables[get_table(k)].displacement(k);
    }

    inline iterator begin()              { return tables[0].begin(); }
    inline iterator end()                { return tables[0].end(); }
    inline const_iterator begin()  const { return tables[0].begin(); }
    inline const_iterator end()    const { return tables[0].end(); }
    inline const_iterator cbegin() const { return tables[0].cbegin(); }
    inline const_iterator cend()   const { return tables[0].cend(); }

private:
    inline size_t get_table(key_type k) const
    {
        return hasher(k) & bits;
    }

public:
    inline static void print_init_header(otm::output_type& out)
    {
        out << otm::width(10) << "f_cap";
    }

    inline void print_init_data(otm::output_type& out)
    {
        size_t cap = 0;
        for (size_t i = 0; i < tl; ++i)
            cap += tables[i].get_capacity();

        out << otm::width(10) << cap;
    }
};


};
