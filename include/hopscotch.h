#pragma once

#include "hopscotch_map.h"
//#include "extern/hopscotch-map/src/hopscotch_map.h"
#include "config.h"
#include "bucket.h"
#include "iterator_base.h"


template<class K, class D, class HF = std::hash<K>, class Conf = Config<> >
class Hopscotch
{
private:
    using Table_t    = tsl::hopscotch_map<K,D,HF>;
    using IterNative_t = typename Table_t::iterator;
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
    using Iterator = typename Table_t::iterator;

    Hopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint)
    { }

    Hopscotch(const Hopscotch&) = delete;
    Hopscotch& operator=(const Hopscotch&) = delete;

    Hopscotch(Hopscotch&&) = default;
    Hopscotch& operator=(Hopscotch&&) = default;

    std::pair<Iterator,bool> insert(Key k, Data d)
    {
        return insert(std::make_pair(k,d));
    }
    std::pair<Iterator,bool> insert(std::pair<Key,Data> t)
    {
        // really stupid shit
        auto ret = table.insert(t);

        return table.insert(t);
    }
    Iterator find  (const Key k)
    {
        return table.find(k);
        //return std::make_pair((it!=end(table)), it->second);
    }
    bool remove(Key k)
    {
        return table.erase(k);
    }

private:
    Table_t table;

public:
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
    static constexpr size_t nh_size = (Conf::NeighborSize <= 62) ?
                                          Conf::NeighborSize : 62;
    using Table_t = tsl::hopscotch_map<K,D,HF, std::equal_to<K>, std::allocator<std::pair<K,D> >,
                                       nh_size, typename Conf::GrowRatio>;
    using IterNative_t = typename Table_t::iterator;

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
    using Pair_t = std::pair<Key,Data>;

    //using Iterator = typename Table_t::iterator;
    class Iterator : public IterNative_t
    {
    public:
        using difference_type   = typename IterNative_t::difference_type;
        using value_type = Pair_t;
        using reference  = Pair_t&;
        using pointer    = Pair_t*;
        using iterator_category = typename IterNative_t::iterator_category;

        Iterator(const IterNative_t& rhs) : IterNative_t(rhs) { };

        reference operator*() const
        {
            return const_cast<reference>(IterNative_t::operator*());
        }
        pointer operator->() const
        {
            return const_cast<pointer>(IterNative_t::operator->());
        }
    };

    SpaceHopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint)
    { }

    SpaceHopscotch(const SpaceHopscotch&) = delete;
    SpaceHopscotch& operator=(const SpaceHopscotch&) = delete;

    SpaceHopscotch(SpaceHopscotch&&) = default;
    SpaceHopscotch& operator=(SpaceHopscotch&&) = default;

    inline std::pair<Iterator, bool> insert(const Key k, const Data d)
    { return insert(std::make_pair(k,d)); }
    inline std::pair<Iterator, bool> insert(std::pair<Key,Data> t)
    { return table.insert(t); }

    inline Iterator find(const Key k)
    { return table.find(k); }
    inline bool remove(const Key k)
    { return table.erase(k); }

private:
    Table_t table;

public:
    inline static void print_init_header(std::ostream& out)
    {
        out.width(6); out << "nghb";
        out.width(6); out << "grat";
        out.width(9); out << "f_cap";
        out << std::flush;
    }
    inline void print_init_data  (std::ostream& out)
    {
        out.width(6); out << nh_size;
        out.width(6); out << Conf::GrowRatio_d;
        out.width(9); out << table.bucket_count();
        out << std::flush;
    }

};
