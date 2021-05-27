#pragma once

/*******************************************************************************
 * include/prob_quadratic.h
 *
 * prob_quadratic and prob_quadratic_inplace implement quadratic probing
 * adhering to the size constraint.  The inplace variant uses memory
 * overcommiting to resize the table without full table reallocations,
 * that would temporarily violate the memory constraint.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.hpp"
#include "prob_base.h"

namespace dysect
{


/* Using Classic Quadratic Probing to Fill a Table Densely ********************/
    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_quadratic : public prob_traits<prob_quadratic<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_quadratic<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;
    public:
        using size_type      = typename base_type::size_type;
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

        prob_quadratic(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
            : base_type(std::max<size_type>(cap, 500), size_constraint)
        {
        }
        prob_quadratic(const prob_quadratic&) = delete;
        prob_quadratic& operator=(const prob_quadratic&) = delete;
        prob_quadratic(prob_quadratic&& rhs)  = default;
        prob_quadratic& operator=(prob_quadratic&& ) = default;

        std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
        std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>&t);
        iterator       find(const key_type& k);
        const_iterator find(const key_type& k) const;
        size_type      erase(const key_type& k);
        int            displacement(const key_type& k) const;
    private:
        using base_type::alpha;
        using base_type::n;
        using base_type::capacity;
        using base_type::table;
        using base_type::h;
        using base_type::inc_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

        static constexpr size_type bitmask = (1ull<<32)-1;


        // Access Functions ****************************************************
        inline size_type index(size_type i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_type mod(size_type i)   const
        {
            auto tmp = i;
            while(tmp >= capacity) tmp -= capacity;
            return tmp;
        }
        inline size_type next(size_t ind, size_t i) const
        {
            return (i < capacity) ? mod(ind +2*i +1) : mod(ind+1);
        }

        // Growing *************************************************************
        inline void grow()
        {
            auto ntable = this_type(n, alpha);

            for (size_type i = 0; i < capacity; ++i)
            {
                auto temp = table[i];
                if (temp.first)
                {
                    ntable.insert(temp);
                }
            }
            (*this) = std::move(ntable);
        }
        // Specialized Deletion stuff ******************************************
    };



    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_quadratic<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_quadratic<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };


    template <class K, class D, class HF, class C>
    inline std::pair<typename prob_quadratic<K,D,HF,C>::iterator, bool>
    prob_quadratic<K,D,HF,C>::insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    template <class K, class D, class HF, class C>
    inline std::pair<typename prob_quadratic<K,D,HF,C>::iterator, bool>
    prob_quadratic<K,D,HF,C>::insert(const std::pair<key_type, mapped_type>&t)
    {
        size_t ind = h(t.first);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if ( temp.first == 0)
            {
                table[ind] = t;
                inc_n();
                return std::make_pair(make_iterator(&table[ind]), true);
            }
            if (temp.first == t.first)
            {
                return std::make_pair(make_iterator(&table[ind]), false);
            }
            ind = next(ind, i);
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::iterator
    prob_quadratic<K,D,HF,C>::find(const key_type& k)
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if (temp.first == k)
            {
                return make_iterator(&table[ind]);
            }
            if ( temp.first == 0)
            {
                return base_type::end();
            }
            ind = next(ind, i);
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::const_iterator
    prob_quadratic<K,D,HF,C>::find(const key_type& k) const
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if (temp.first == k)
            {
                return make_citerator(&table[ind]);
            }
            if ( temp.first == 0)
            {
                return base_type::cend();
            }
            ind = next(ind, i);

        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::size_type
    prob_quadratic<K,D,HF,C>::erase([[maybe_unused]] const key_type& k)
    {
        // TODO
        return 0;
    }

    template <class K, class D, class HF, class C>
    inline int
    prob_quadratic<K,D,HF,C>::displacement(const key_type& k) const
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            if (table[ind].first == k)
            {
                return i;
            }
            ind = next(ind, i);

        }
    }




















// *****************************************************************************
// Same as Above, but Growing Using in Place Migration *************************
// *****************************************************************************

    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_quadratic_inplace : public prob_traits<prob_quadratic_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_quadratic_inplace<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;

    public:
        using size_type      = typename base_type::size_type;
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;
    private:
        using value_intern = std::pair<key_type, mapped_type>;

        static constexpr size_type max_size = 16ull << 30;

    public:
        prob_quadratic_inplace(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
            : base_type(0, size_constraint), max_buffer_size(0)
        {
            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            if (temp) table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            thresh   = (cap) ? cap*beta  : 2048*beta;

            std::fill(table.get(), table.get()+capacity, value_intern());
        }
        prob_quadratic_inplace(const prob_quadratic_inplace&) = delete;
        prob_quadratic_inplace& operator=(const prob_quadratic_inplace&) = delete;
        prob_quadratic_inplace(prob_quadratic_inplace&& rhs)  = default;
        prob_quadratic_inplace& operator=(prob_quadratic_inplace&& ) = default;

        std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
        std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>&t);
        iterator       find(const key_type& k);
        const_iterator find(const key_type& k) const;
        size_type      erase(const key_type& k);
        int            displacement(const key_type& k) const;
    private:
        using base_type::alpha;
        using base_type::beta;
        using base_type::n;
        using base_type::capacity;
        using base_type::thresh;
        using base_type::table;
        using base_type::inc_n;
        using base_type::dec_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

        size_type max_buffer_size;

        static constexpr size_type bitmask = (1ull<<32)-1;

        using base_type::h;

        // Access Functions ****************************************************
        inline size_type index(size_type i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_type mod(size_type i)   const
        {
            auto tmp = i;
            while(tmp >= capacity) tmp -= capacity;
            return tmp;
        }
        inline size_type next(size_t ind, size_t i) const
        {
            return (i < capacity) ? mod(ind +2*i +1) : mod(ind+1);
        }

    private:
        // Growing *************************************************************
        inline void grow()
        {
            size_type ncap    = n*alpha;
            size_type nthresh = n*beta;
            std::fill(table.get()+capacity, table.get()+ncap, value_intern());

            auto ocap = capacity;
            capacity = ncap;
            thresh   = nthresh;

            std::vector<value_intern> buffer;

            for (int i = ocap -1; i >= 0; --i)
            {
                auto temp = table[i];

                if (temp.first)
                {
                    table[i] = value_intern();

                    if (!intern_migration_insert(i, temp))
                        buffer.push_back(temp);
                }
            }

            max_buffer_size = std::max(max_buffer_size, buffer.size());
            for (auto it = buffer.begin(); it != buffer.end(); it++)
            {
                dec_n();
                insert(*it);
            }
        }

        inline bool intern_migration_insert(uint low, const value_intern& e)
        {
            auto ind = h(e.first);

            // x+(i+1)² = x+i²+2i+1²
            for (size_t i=0; ind>low; ++i)
            {
                if (table[ind].first == 0)
                {
                    table[ind] = e;
                    return true;
                }
                ind = next(ind, i);
            }
            return false;
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_quadratic_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_quadratic_inplace<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };


    template <class K, class D, class HF, class C>
    inline std::pair<typename prob_quadratic_inplace<K,D,HF,C>::iterator, bool>
    prob_quadratic_inplace<K,D,HF,C>::insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    template <class K, class D, class HF, class C>
    inline std::pair<typename prob_quadratic_inplace<K,D,HF,C>::iterator, bool>
    prob_quadratic_inplace<K,D,HF,C>::insert(const std::pair<key_type, mapped_type>&t)
    {
        size_t ind = h(t.first);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if ( temp.first == 0)
            {
                table[ind] = t;
                inc_n();
                return std::make_pair(make_iterator(&table[ind]), true);
            }
            if (temp.first == t.first)
            {
                return std::make_pair(make_iterator(&table[ind]), false);
            }
            ind = next(ind, i);
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic_inplace<K,D,HF,C>::iterator
    prob_quadratic_inplace<K,D,HF,C>::find(const key_type& k)
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if (temp.first == k)
            {
                return make_iterator(&table[ind]);
            }
            if ( temp.first == 0)
            {
                return base_type::end();
            }
            ind = next(ind, i);
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic_inplace<K,D,HF,C>::const_iterator
    prob_quadratic_inplace<K,D,HF,C>::find(const key_type& k) const
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            auto temp  = table[ind];

            if (temp.first == k)
            {
                return make_citerator(&table[ind]);
            }
            if ( temp.first == 0)
            {
                return base_type::cend();
            }
            ind = next(ind, i);
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic_inplace<K,D,HF,C>::size_type
    prob_quadratic_inplace<K,D,HF,C>::erase(const key_type& k)
    {
        // TODO
        return 0;
    }

    template <class K, class D, class HF, class C>
    inline int
    prob_quadratic_inplace<K,D,HF,C>::displacement(const key_type& k) const
    {
        size_t ind = h(k);

        // x+(i+1)² = x+i²+2i+1²
        for (size_t i=0; ; ++i)
        {
            if (table[ind].first == k)
            {
                return i;
            }
            ind = next(ind, i);
        }
    }

} // namespace dysect
