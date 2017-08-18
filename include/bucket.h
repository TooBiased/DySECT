#pragma once

/*******************************************************************************
 * include/bucket.h
 *
 * bucket implementation for variants of bucket cuckoo hashing.
 * This implementation ensures that all contained elements are stored in the
 * beginning of the bucket. This allows some performance tricks.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tuple>



template<class K, class D, size_t BS = 4>
class Bucket
{
public:
    using Key = K;
    using Data = D;
    using Pair_t = std::pair<Key,Data>;
    using FRet = std::pair<bool, Data>;

    Bucket() { for (size_t i = 0; i < BS; ++i) elements[i] = std::make_pair(Key(), Data());}
    Bucket(const Bucket& rhs) = default;
    Bucket& operator=(const Bucket& rhs) = default;

    bool   insert(const Key& k, const Data& d);
    bool   insert(const Pair_t& t);
    FRet   find  (const Key& k);
    bool   remove(const Key& k);
    FRet   pop   (const Key& k);

    int    probe (const Key& k);

    bool   space ();
    Pair_t get(const size_t i);
    Pair_t replace(const size_t i, const Pair_t& t);


    Pair_t* insertPtr(const Pair_t& t);
    const Pair_t* findPtr(const Key& k) const;
    Pair_t* findPtr(const Key& k);
    std::pair<int, Pair_t*> probePtr(const Key& k);

    Pair_t elements[BS];
};


template<class K, class D, size_t BS>
inline bool Bucket<K,D,BS>::insert(const Key& k, const Data& d)
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
inline bool Bucket<K,D,BS>::insert(const Pair_t& t)
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
inline typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::find(const Key& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first )      return std::make_pair(false, Data());
        if ( elements[i].first == k ) return std::make_pair(true, elements[i].second);
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
inline bool Bucket<K,D,BS>::remove(const Key& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            size_t j = BS-1;
            for ( ; !elements[j].first; --j ) { }
            elements[i] = elements[j];
            elements[j] = std::make_pair(Key(), Data());
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
inline typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::pop(const Key& k)
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
inline int Bucket<K,D,BS>::probe(const Key& k)
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
inline std::pair<K, D> Bucket<K,D,BS>::replace(size_t i, const Pair_t& newE)
{
    auto temp = elements[i];
    elements[i] = newE;
    return temp;
}


template<class K, class D, size_t BS>
inline std::pair<K,D>* Bucket<K,D,BS>::insertPtr(const Pair_t& t)
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
inline std::pair<K,D>* Bucket<K,D,BS>::findPtr(const Key& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        //if (!elements[i].first )      return nullptr;
        if ( elements[i].first == k ) return &elements[i];
    }
    return nullptr;
}

template<class K, class D, size_t BS>
inline const std::pair<K,D>* Bucket<K,D,BS>::findPtr(const Key& k) const
{
    for (size_t i = 0; i < BS; ++i)
    {
        //if (!elements[i].first )      return nullptr;
        if ( elements[i].first == k ) return &elements[i];
    }
    return nullptr;
}

template<class K, class D, size_t BS>
inline std::pair<int, std::pair<K,D>*> Bucket<K,D,BS>::probePtr(const Key& k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first)      return std::make_pair(BS - i, &elements[i]);
        if ( elements[i].first == k) return std::make_pair(-1    , &elements[i]);
    }
    return std::make_pair(0, nullptr);
}
