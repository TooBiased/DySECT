#pragma once

/*******************************************************************************
 * include/bucket.h
 *
 * bucket implementation for variants of bucket cuckoo hashing.
 * This implementation ensures that all contained elements are stored in the
 * beginning of the bucket. This allows some performance optimizations.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tuple>

namespace dysect
{

    template<class K, class D, size_t BS = 4>
    class bucket
    {
    public:
        using key_type    = K;
        using mapped_type = D;
    private:
        using value_intern     = std::pair<key_type,mapped_type>;
    public:
        using find_return_type = std::pair<bool, mapped_type>;

        bucket()
        {
            for (size_t i = 0; i < BS; ++i)
                elements[i] = value_intern();
        }
        bucket(const bucket& rhs) = default;
        bucket& operator=(const bucket& rhs) = default;

        bool   insert(const key_type& k, const mapped_type& d);
        bool   insert(const value_intern& t);
        find_return_type   find  (const key_type& k);
        bool   remove(const key_type& k);
        find_return_type   pop   (const key_type& k);

        int    probe (const key_type& k);
        int    displacement(const key_type& k) const;

        bool   space ();
        value_intern get(const size_t i);
        value_intern replace(const size_t i, const value_intern& t);


        value_intern* insert_ptr(const value_intern& t);
        const value_intern* find_ptr(const key_type& k) const;
        value_intern* find_ptr(const key_type& k);
        std::pair<int, value_intern*> probe_ptr(const key_type& k);

        value_intern elements[BS];
    };


    template<class K, class D, size_t BS>
    inline bool bucket<K,D,BS>::insert(const key_type& k, const mapped_type& d)
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
    inline bool bucket<K,D,BS>::insert(const value_intern& t)
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
    inline typename bucket<K,D,BS>::find_return_type bucket<K,D,BS>::find(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if (!elements[i].first )      return std::make_pair(false, mapped_type());
            if ( elements[i].first == k ) return std::make_pair(true, elements[i].second);
        }
        return std::make_pair(false, mapped_type());
    }

    template<class K, class D, size_t BS>
    inline bool bucket<K,D,BS>::remove(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if (elements[i].first == k)
            {
                size_t j = BS-1;
                for ( ; !elements[j].first; --j ) { }
                elements[i] = elements[j];
                elements[j] = std::make_pair(key_type(), mapped_type());
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
    inline typename bucket<K,D,BS>::find_return_type bucket<K,D,BS>::pop(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if (elements[i].first == k)
            {
                mapped_type d = elements[i].second;
                for (size_t j = i+1; j < BS; ++j)
                {
                    if (elements[j].first) { elements[i] = elements[j]; i = j; }
                    else break;
                }
                elements[i] = std::make_pair(key_type(), mapped_type());
                return std::make_pair(true, d);
            }
            else if (! elements[i].first)
            {
                break;
            }
        }
        return std::make_pair(false, mapped_type());
    }

    template<class K, class D, size_t BS>
    inline int bucket<K,D,BS>::probe(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if (!elements[i].first)      return BS - i;
            if ( elements[i].first == k) return -1;

        }
        return 0;
    }

    template<class K, class D, size_t BS>
    inline int bucket<K,D,BS>::displacement(const key_type& k) const
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if ( elements[i].first == k ) return i;
        }
        return BS;
    }

    template<class K, class D, size_t BS>
    inline bool bucket<K,D,BS>::space()
    {
        return !elements[BS-1].first;
    }

    template<class K, class D, size_t BS>
    inline std::pair<K, D> bucket<K,D,BS>::get(size_t i)
    {
        return elements[i];
    }

    template<class K, class D, size_t BS>
    inline std::pair<K, D> bucket<K,D,BS>::replace(size_t i, const value_intern& newE)
    {
        auto temp = elements[i];
        elements[i] = newE;
        return temp;
    }


    template<class K, class D, size_t BS>
    inline std::pair<K,D>* bucket<K,D,BS>::insert_ptr(const value_intern& t)
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
    inline std::pair<K,D>* bucket<K,D,BS>::find_ptr(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            //if (!elements[i].first )      return nullptr;
            if ( elements[i].first == k ) return &elements[i];
        }
        return nullptr;
    }

    template<class K, class D, size_t BS>
    inline const std::pair<K,D>* bucket<K,D,BS>::find_ptr(const key_type& k) const
    {
        for (size_t i = 0; i < BS; ++i)
        {
            //if (!elements[i].first )      return nullptr;
            if ( elements[i].first == k ) return &elements[i];
        }
        return nullptr;
    }

    template<class K, class D, size_t BS>
    inline std::pair<int, std::pair<K,D>*> bucket<K,D,BS>::probe_ptr(const key_type& k)
    {
        for (size_t i = 0; i < BS; ++i)
        {
            if (!elements[i].first)      return std::make_pair(BS - i, &elements[i]);
            if ( elements[i].first == k) return std::make_pair(-1    , &elements[i]);
        }
        return std::make_pair(0, nullptr);
    }

} // namespace dysect
