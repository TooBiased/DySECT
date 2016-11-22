#pragma once

#include "hopscotch_map.h"
//#include "extern/hopscotch-map/src/hopscotch_map.h"
#include "cuckoo_base.h"


template<class K, class D, class HF = std::hash<K>, class Config = CuckooConfig<> >
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
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint), hcounter(0)
    { }

    Hopscotch(const Hopscotch&) = delete;
    Hopscotch& operator=(const Hopscotch&) = delete;

    Hopscotch(Hopscotch&&) = default;
    Hopscotch& operator=(Hopscotch&&) = default;

    bool insert(Key k, Data d) { return insert(std::make_pair(k,d)); } //make bool
    bool insert(std::pair<Key,Data> t)
    {
        return table.insert(t).second;
    }
    FRet find  (Key k) const
    {
        auto it = table.find(k);
        return std::make_pair((it!=end(table)), it->second);
    } //correct
    bool remove(Key k)
    {
        return false;
    }

private:
    Table_t table;

public:
    /*** for symmetry with my implementation **********************************/
    using HistCount_t = typename Config::HistCount_t;
    using Bucket_t    = Bucket<K,D,1>;
    HistCount_t hcounter;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    std::pair<size_t, Bucket_t*> getTable(size_t) { return std::make_pair(0ul, nullptr); }
    void clearHist() { }
};


template<class K, class D, class HF = std::hash<K>, class Config = CuckooConfig<> >
class SpaceHopscotch
{
private:
    using Table_t = tsl::hopscotch_map<K,D,HF, std::equal_to<K>, std::allocator<std::pair<K,D> >, 32, std::ratio<5,4> >;
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

    SpaceHopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint), hcounter(0)
    { }

    SpaceHopscotch(const SpaceHopscotch&) = delete;
    SpaceHopscotch& operator=(const SpaceHopscotch&) = delete;

    SpaceHopscotch(SpaceHopscotch&&) = default;
    SpaceHopscotch& operator=(SpaceHopscotch&&) = default;

    bool insert(Key k, Data d) { return insert(std::make_pair(k,d)); } //make bool
    bool insert(std::pair<Key,Data> t)
    {
        return table.insert(t).second;
    }
    FRet find  (Key k) const
    {
        auto it = table.find(k);
        return std::make_pair((it!=table.end()), it->second);
    } //correct
    bool remove(Key k)
    {
        return false;
    }

private:
    Table_t table;

public:
    /*** for symmetry with my implementation **********************************/
    using HistCount_t = typename Config::HistCount_t;
    using Bucket_t    = Bucket<K,D,1>;
    HistCount_t hcounter;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    std::pair<size_t, Bucket_t*> getTable(size_t) { return std::make_pair(0ul, nullptr); }
    void clearHist() { }
};
