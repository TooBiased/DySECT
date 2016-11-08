#pragma once

#include <tuple>

template<class K, class D, size_t BS = 4>
class Bucket
{
public:
    using Key = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    Bucket() { for (size_t i = 0; i < BS; ++i) elements[i] = std::make_pair(Key(), Data());}
    Bucket(const Bucket& rhs) = default;
    Bucket& operator=(const Bucket& rhs) = default;

    bool   insert(Key k, Data d);
    bool   insert(std::pair<Key, Data> t);
    FRet   find  (Key k);
    bool   remove(Key k);
    FRet   pop   (Key k);

    int    probe (Key k);

    bool   space ();
    std::pair<Key, Data> get(size_t i);
    std::pair<Key, Data> replace(size_t i, std::pair<Key, Data> t);

    std::pair<Key, Data> elements[BS];
};


template<class K, class D, size_t BS>
inline bool Bucket<K,D,BS>::insert(Key k, Data d)
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
inline bool Bucket<K,D,BS>::insert(std::pair<Key,Data> t)
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
inline typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::find(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first )      return std::make_pair(false, Data());
        if ( elements[i].first == k ) return std::make_pair(true, elements[i].second);
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
inline bool Bucket<K,D,BS>::remove(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            for (size_t j = i+1; j < BS; ++j)
            {
                if (elements[j].first) { elements[i] = elements[j]; i = j; }
                else break;
            }
            elements[i] = std::make_pair(Key(), Data());
            return true;
        }
        else if (! elements[i].first)
        {
            break;
        }
    }
    return false;
}

template<class K, class D, size_t BS>
inline typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::pop(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            Data d = elements[i].second;
            for (size_t j = i+1; j < BS; ++j)
            {
                if (elements[j].first) { elements[i] = elements[j]; i = j; }
                else break;
            }
            elements[i] = std::make_pair(Key(), Data());
            return std::make_pair(true, d);
        }
        else if (! elements[i].first)
        {
            break;
        }
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
inline int Bucket<K,D,BS>::probe(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first)      return BS - i;
        if ( elements[i].first == k) return -1;

    }
    return 0;
}

template<class K, class D, size_t BS>
inline bool Bucket<K,D,BS>::space()
{
    return !elements[BS-1].first;
}

template<class K, class D, size_t BS>
inline std::pair<K, D> Bucket<K,D,BS>::get(size_t i)
{
    return elements[i];
}

template<class K, class D, size_t BS>
inline std::pair<K, D> Bucket<K,D,BS>::replace(size_t i, std::pair<K, D> newE)
{
    auto temp = elements[i];
    elements[i] = newE;
    return temp;
}
