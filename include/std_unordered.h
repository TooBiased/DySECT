#pragma once

#include <unordered_map>

template <class K, class D, class HF = std::hash<K>,
          class Conf = Config<> >
class STDProb
{
private:
    using Table_t  = std::unordered_map<K,D,HF>;
    using iterator = typename Table_t::iterator;

    Table_t table;

public:
    using key_type    = K;
    using mapped_type = D;

    STDProb(size_t cap, double = 1.1, size_t = 0) : table(cap) {}

    inline std::pair<iterator, bool> insert(const std::pair<key_type,mapped_type> t)
    {
        return table.insert(t);
    }

    inline std::pair<iterator, bool> insert(const key_type k, const mapped_type d)
    {
        return table.insert(std::make_pair(k, d));
    }

    inline iterator find(const key_type k)
    {
        return table.find(k);
    }

    inline size_t erase(const key_type k)
    {
        return table.erase(k);
    }

    inline iterator       begin()  const { return table.begin(); }
    inline iterator       end()    const { return table.end();   }
    inline const_iterator cbegin() const { return table.cbegin(); }
    inline const_iterator cend()   const { return table.cend();   }

    inline static void print_init_header(std::ostream& )
    { }
    inline void print_init_data(std::ostream& )
    { }
};
