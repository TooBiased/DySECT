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
    // template<class key_type,
    //          class T,
    //          class Hash = std::hash<key_type>,
    //          class key_typeEqual = std::equal_to<key_type>,
    //          class Allocator = std::allocator<std::pair<key_type, T>>,
    //          unsigned int NeighborhoodSize = 62,
    //          class GrowthFactor = std::ratio<2, 1>>
    // class hopscotch_map;


    using IterNative_t = typename Table_t::iterator;
public:
    using key_type  = K;
    using mapped_type = D;
    using value_intern = std::pair<key_type,mapped_type>;
    class iterator : public IterNative_t
    {
    public:
        using difference_type   = typename IterNative_t::difference_type;
        using value_type = std::pair<const key_type,mapped_type>;
        using reference  = value_type&;
        using pointer    = value_type*;
        using iterator_category = typename IterNative_t::iterator_category;

        iterator(const IterNative_t& rhs) : IterNative_t(rhs) { };

        reference operator*() const
        {
            return reinterpret_cast<reference>(
                       const_cast<value_intern&>(IterNative_t::operator*()));
        }
        pointer operator->() const
        {
            return reinterpret_cast<pointer>(
                       const_cast<value_intern*>(IterNative_t::operator->()));
        }
    };
    using const_iterator = typename Table_t::const_iterator;

    Hopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint)
    { }

    Hopscotch(const Hopscotch&) = delete;
    Hopscotch& operator=(const Hopscotch&) = delete;

    Hopscotch(Hopscotch&&) = default;
    Hopscotch& operator=(Hopscotch&&) = default;

    inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
    { return insert(std::make_pair(k,d)); }
    inline std::pair<iterator, bool> insert(const std::pair<key_type,mapped_type>& t)
    { return table.insert(t); }

    inline iterator find(const key_type& k)
    { return table.find(k); }
    inline const_iterator find(const key_type& k) const
    { return table.find(k); }
    inline size_t erase(const key_type& k)
    { return table.erase(k); }

    inline iterator begin()  const { return table.begin();  }
    inline iterator end()    const { return table.end();    }
    inline iterator cbegin() const { return table.cbegin(); }
    inline iterator cend()   const { return table.cend();   }

    inline mapped_type& at(const key_type& k)             { return table.at(k); }
    inline const mapped_type& at(const key_type& k) const { return table.at(k); }
    inline mapped_type& operator[](const key_type& k)     { return table[k]; }
    inline size_t count(const key_type& k)          const { return table.count(k); }

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
    using mapped_type = D;
    using value_intern = std::pair<Key,mapped_type>;

    //using iterator = typename Table_t::iterator;
    class iterator : public IterNative_t
    {
    public:
        using difference_type   = typename IterNative_t::difference_type;
        using value_type = std::pair<const Key,mapped_type>;
        using reference  = value_type&;
        using pointer    = value_type*;
        using iterator_category = typename IterNative_t::iterator_category;

        iterator(const IterNative_t& rhs) : IterNative_t(rhs) { };

        reference operator*() const
        {
            return reinterpret_cast<reference>(
                       const_cast<value_intern&>(IterNative_t::operator*()));
        }
        pointer operator->() const
        {
            return reinterpret_cast<pointer>(
                       const_cast<value_intern*>(IterNative_t::operator->()));
        }
    };
    using const_iterator = typename Table_t::const_iterator;

    SpaceHopscotch(size_t cap = 0, double size_constraint = 1.1,
              size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
        : table(cap*size_constraint)
    { }

    SpaceHopscotch(const SpaceHopscotch&) = delete;
    SpaceHopscotch& operator=(const SpaceHopscotch&) = delete;

    SpaceHopscotch(SpaceHopscotch&&) = default;
    SpaceHopscotch& operator=(SpaceHopscotch&&) = default;

    inline std::pair<iterator, bool> insert(const Key& k, const mapped_type& d)
    { return insert(std::make_pair(k,d)); }
    inline std::pair<iterator, bool> insert(const std::pair<Key,mapped_type>& t)
    { return table.insert(t); }

    inline iterator find(const Key& k)
    { return table.find(k); }
    inline const_iterator find(const Key& k) const
    { return table.find(k); }
    inline size_t erase(const Key& k)
    { return table.erase(k); }

    inline iterator begin()        { return table.begin();  }
    inline iterator end()          { return table.end();    }
    inline iterator cbegin() const { return table.cbegin(); }
    inline iterator cend()   const { return table.cend();   }

    inline mapped_type& at(const Key& k)             { return table.at(k);    }
    inline const mapped_type& at(const Key& k) const { return table.at(k);    }
    inline mapped_type& operator[](const Key& k)     { return table[k];       }
    inline size_t count(const Key& k)   const { return table.count(k); }

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
