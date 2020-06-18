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
            factor = double(capacity)/double(1ull << 32);
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

        double factor;

        static constexpr size_type bitmask = (1ull<<32)-1;


        // Access Functions ****************************************************
        inline size_type index(size_type i) const
        { return (bitmask & i)*factor; }
        inline size_type mod(size_type i)   const
        { return (i < capacity) ? i : i-capacity; }

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
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if ( temp.first == 0)
            {
                table[ti] = t;
                inc_n();
                return std::make_pair(make_iterator(&table[ti]), true);
            }
            if (temp.first == t.first)
            {
                return std::make_pair(make_iterator(&table[ti]), false);
            }
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::iterator
    prob_quadratic<K,D,HF,C>::find(const key_type& k)
    {
        size_t ind = h(k);
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if (temp.first == k)
            {
                return make_iterator(&table[ti]);
            }
            if ( temp.first == 0)
            {
                return base_type::end();
            }
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::const_iterator
    prob_quadratic<K,D,HF,C>::find(const key_type& k) const
    {
        size_t ind = h(k);
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if (temp.first == k)
            {
                return make_citerator(&table[ti]);
            }
            if ( temp.first == 0)
            {
                return base_type::cend();
            }
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic<K,D,HF,C>::size_type
    prob_quadratic<K,D,HF,C>::erase(const key_type& k)
    {
        // TODO
        return 0;
    }

    template <class K, class D, class HF, class C>
    inline int
    prob_quadratic<K,D,HF,C>::displacement(const key_type& k) const
    {
        size_t ind = h(k);
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);

            if (table[ti].first == k)
            {
                return i;
            }
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
            : base_type(0, size_constraint), bla(0)
        {
            // factor = double(capacity-300)/double(1ull << 32);

            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            if (temp) table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            thresh   = (cap) ? cap*beta  : 2048*beta;
            factor = double(capacity)/double(1ull << 32);

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
        using base_type::make_iterator;
        using base_type::make_citerator;

        double    factor;
        size_type bla;

        static constexpr size_type bitmask = (1ull<<32)-1;

        using base_type::h;

        // Access Functions ****************************************************
        inline size_type index(size_type i) const
        { return (bitmask & i)*factor; }
        inline size_type mod(size_type i)   const
        { return (i < capacity) ? i : i-capacity; }

    private:
        // Growing *************************************************************
        inline void grow()
        {
            // TODO
        }

    public:
        inline static void print_init_header(std::ostream& out)
        {
            out.width(6); out  << "busize" << " ";
            base_type::print_init_header(out);
        }

        inline        void print_init_data  (std::ostream& out)
        {
            out.width(6);  out << bla << " ";
            base_type::print_init_data(out);
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
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if ( temp.first == 0)
            {
                table[ti] = t;
                inc_n();
                return std::make_pair(make_iterator(&table[ti]), true);
            }
            if (temp.first == t.first)
            {
                return std::make_pair(make_iterator(&table[ti]), false);
            }
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic_inplace<K,D,HF,C>::iterator
    prob_quadratic_inplace<K,D,HF,C>::find(const key_type& k)
    {
        size_t ind = h(k);
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if (temp.first == k)
            {
                return make_iterator(&table[ti]);
            }
            if ( temp.first == 0)
            {
                return base_type::end();
            }
        }
    }

    template <class K, class D, class HF, class C>
    inline typename prob_quadratic_inplace<K,D,HF,C>::const_iterator
    prob_quadratic_inplace<K,D,HF,C>::find(const key_type& k) const
    {
        size_t ind = h(k);
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);
            auto temp  = table[ti];

            if (temp.first == k)
            {
                return make_citerator(&table[ti]);
            }
            if ( temp.first == 0)
            {
                return base_type::cend();
            }
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
        for (size_t i=0; ; ++i)
        {
            auto ti    = mod(i*i + ind);

            if (table[ti].first == k)
            {
                return i;
            }
        }
    }

} // namespace dysect
