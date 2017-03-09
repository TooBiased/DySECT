#pragma once

#include <unordered_map>

template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class STDProb
{
private:
    using Table_t  = std::unordered_map<K,D,HF>;
    using Iterator = typename Table_t::iterator;

    Table_t table;

public:
    using Key = K;
    using Data = D;

    FastLinProb(size_t cap, double = 1.1, size_t = 0) : table(cap) {}

    inline std::pair<Iterator, bool> insert(const std::pair<Key,Data> t)
    {
        return table.insert(t);
    }

    inline std::pair<Iterator, bool> insert(const Key k, const Data d)
    {
        return table.insert(std::make_pair(k, d));
    }

    inline Iterator find(const Key k)
    {
        return table.find(k);
    }

    inline bool remove(const Key k);

    inline static void print_init_header(std::ostream& )
    { }
    inline void print_init_data(std::ostream& )
    { }

    constexpr Iterator end() { return table.end(); }
};
