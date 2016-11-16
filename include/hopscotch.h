#pragma once

#include "hopscotch_map.h"

template<class K, class D, class HF = std::hash<K>, class Config = size_t>
class Hopscotch
{
private:
    using Table_t = tsl::hopscotch_map<K,D,HF>;
    // template<class Key,
    //          class T,
    //          class Hash = std::hash<Key>,
    //          class KeyEqual = std::equal_to<Key>,
    //          class Allocator = std::allocator<std::pair<Key, T>>,
    //          unsigned int NeighborhoodSize = 62,
    //          class GrowthFactor = std::ratio<2, 1>>
    // class hopscotch_map;

public:
    using Key  = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    Hopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/, size_t /*seed*/ = 0)
        : table(cap*size_constraint)
    { }

    Hopscotch(const Hopscotch&) = delete;
    Hopscotch& operator=(const Hopscotch&) = delete;

    Hopscotch(Hopscotch&&) = default;
    Hopscotch& operator=(Hopscotch&&) = default;

    bool insert(Key k, Data d) { table[k] = d; return true; } //make bool
    bool insert(std::pair<Key,Data> t) { table.insert(t); return true; };
    FRet find  (Key k) const   { std::make_pair(true, table[k]); } //correct
    bool remove(Key k);

private:
    Table_t table;
};
