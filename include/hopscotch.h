#pragma once

#include "hopscotch_map.h"
//#include "extern/hopscotch-map/src/hopscotch_map.h"
#include "config.h"
#include "bucket.h"


template<class K, class D, class HF = std::hash<K>, class Conf = Config<> >
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
        return erase(k);
    }

private:
    Table_t table;

public:
    /*** for symmetry with my implementation **********************************/

    using HistCount_t = typename Conf::HistCount_t;
    using Bucket_t    = Bucket<K,D,1>;
    HistCount_t hcounter;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    std::pair<size_t, Bucket_t*> getTable(size_t) { return std::make_pair(0ul, nullptr); }
    void clearHist() { }

    /* VISUALISATION **********************************************************/
    static void print_init_header(std::ostream&)
        { }
    static void print_init_data(std::ostream&)
        { }

};




template<class K, class D, class HF = std::hash<K>, class Conf = HopscotchConfig<> >
class SpaceHopscotch
{
private:
    using Table_t = tsl::hopscotch_map<K,D,HF, std::equal_to<K>, std::allocator<std::pair<K,D> >,
                                       Conf::NeighborSize, typename Conf::GrowRatio>;
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
        : table(cap*size_constraint)
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
    }
    bool remove(Key k)
    {
        return erase(k);
    }

private:
    Table_t table;

public:
    /*** for symmetry with my implementation **********************************/
    /*
    using HistCount_t = typename Conf::HistCount_t;
    using Bucket_t    = Bucket<K,D,1>;
    HistCount_t hcounter;

    static constexpr size_t bs = 0;
    static constexpr size_t tl = 0;

    std::pair<size_t, Bucket_t*> getTable(size_t) { return std::make_pair(0ul, nullptr); }
    void clearHist() { }
    */

    /* VISUALISATION **********************************************************/
    inline static void print_init_header(std::ostream& out)
    {
        out.width(6); out << "nghb";
        out.width(6); out << "grat";
        out << std::flush;
    }
    inline static void print_init_data  (std::ostream& out)
    {
        out.width(6); out << Conf::NeighborSize;
        out.width(6); out << Conf::GrowRatio_d;
        out << std::flush;
    }

};
