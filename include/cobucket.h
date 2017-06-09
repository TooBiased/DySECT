#pragma once

#include <tuple>



template<class K, class D, size_t BS = 4>
class co_bucket
{
private:
    using this_type        = co_bucket<K,D,BS>;
public:
    using key_type         = K;
    using mapped_type      = D;
    using value_intern     = std::pair<key_type,mapped_type>;
    using find_return_type = std::pair<bool, mapped_type>;

    co_bucket() { for (size_t i = 0; i < BS; ++i) elements[i] = std::make_pair(key_type(), mapped_type());}
    co_bucket(const co_bucket& rhs) = default;
    co_bucket& operator=(const co_bucket& rhs) = default;

    bool   insert(const key_type& k, const mapped_type& d);
    bool   insert(const value_intern& t);
    find_return_type   find  (const key_type& k);
    bool   remove(const key_type& k);
    find_return_type   pop   (const key_type& k);

    int    probe (const key_type& k);

    bool   space ();
    value_intern get(const size_t i);
    value_intern replace(const size_t i, const value_intern& t);


    value_intern* insertPtr(const value_intern& t);
    const value_intern* findPtr(const key_type& k) const;
    value_intern* findPtr(const key_type& k);
    std::pair<int, value_intern*> probePtr(const key_type& k);

    value_intern elements[BS];
};


template<class K, class D, size_t BS>
inline bool co_bucket<K,D,BS>::insert(const key_type& k, const mapped_type& d)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first) continue;
        elements[i].first  = k;
        elements[i].second = d;

        return true;
    }

    return false;
}

template<class K, class D, size_t BS>
inline bool co_bucket<K,D,BS>::insert(const value_intern& t)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first) continue;
        elements[i]  = t;
        return true;
    }

    return false;
}

template<class K, class D, size_t BS>
inline typename co_bucket<K,D,BS>::find_return_type co_bucket<K,D,BS>::find(const key_type& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if ( elements[i].first == k ) return std::make_pair(true, elements[i].second);
    }
    return std::make_pair(false, mapped_type());
}

template<class K, class D, size_t BS>
inline bool co_bucket<K,D,BS>::remove(const key_type& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            elements[i] = std::make_pair(key_type(), mapped_type());
            return true;
        }
    }
    return false;
}

template<class K, class D, size_t BS>
inline typename co_bucket<K,D,BS>::find_return_type co_bucket<K,D,BS>::pop(const key_type& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            mapped_type d = elements[i].second;
            elements[i] = std::make_pair(key_type(), mapped_type());
            return std::make_pair(true, d);
        }
    }
    return std::make_pair(false, mapped_type());
}

template<class K, class D, size_t BS>
inline int co_bucket<K,D,BS>::probe(const key_type& k)
{
    size_t count = 0;
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first) ++count;
        if ( elements[i].first == k) return -1;

    }
    return count;
}

template<class K, class D, size_t BS>
inline bool co_bucket<K,D,BS>::space()
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first) return true;
    }
    return false;
}

template<class K, class D, size_t BS>
inline std::pair<K, D> co_bucket<K,D,BS>::get(size_t i)
{
    return elements[i];
}

template<class K, class D, size_t BS>
inline std::pair<K, D> co_bucket<K,D,BS>::replace(size_t i, const value_intern& newE)
{
    auto temp = elements[i];
    elements[i] = newE;
    return temp;
}


template<class K, class D, size_t BS>
inline std::pair<K,D>* co_bucket<K,D,BS>::insertPtr(const value_intern& t)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first) continue;

        elements[i]  = t;
        return &elements[i];
    }

    return nullptr;
}

template<class K, class D, size_t BS>
inline std::pair<K,D>* co_bucket<K,D,BS>::findPtr(const key_type& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if ( elements[i].first == k ) return &elements[i];
    }
    return nullptr;
}

template<class K, class D, size_t BS>
inline const std::pair<K,D>* co_bucket<K,D,BS>::findPtr(const key_type& k) const
{
    for (size_t i = 0; i < BS; ++i)
    {
        if ( elements[i].first == k ) return &elements[i];
    }
    return nullptr;
}

template<class K, class D, size_t BS>
inline std::pair<int, std::pair<K,D>*> co_bucket<K,D,BS>::probePtr(const key_type& k)
{
    size_t        count = 0;
    value_intern* tptr  = nullptr;
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first)
        {
            ++count;
            tptr = (tptr) ? tptr : &elements[i];
        }
        if ( elements[i].first == k)
            return std::make_pair(-1, &elements[i]);
    }
    return std::make_pair(count, tptr);
}
