#pragma once

/*******************************************************************************
 * include/prob_simple.h
 *
 * fast_lin_prob implements a fast linear probing hash table, that
 * uses bitmasks for table offsets.  This variant does not respect the
 * size constraint in the same way other tables do.
 *
 * lin_prob and lin_prob_in_place implement linear probing adhering to
 * the size constraint.  The inplace variant uses memory overcommiting
 * to resize the table without full table reallocations, that would
 * temporarily violate the memory constraint.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.hpp"

#include "prob_base.h"

namespace dysect
{


/* Fast linear probing table using powers of 2 and bitmasking *****************/
    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_linear_doubling : public prob_traits<prob_linear_doubling<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_linear_doubling<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;

    public:
        using size_type   = typename base_type::size_type;
        using key_type    = typename prob_traits<this_type>::key_type;
        using mapped_type = typename prob_traits<this_type>::mapped_type;


        prob_linear_doubling(size_type cap = 0, double = 1., size_type = 0)
            : base_type(next_power_of_two(cap) << 1, 1.)
            {
                thresh  = capacity * .6;
                bitmask = capacity - 1;
            }
        prob_linear_doubling(const prob_linear_doubling&) = delete;
        prob_linear_doubling& operator=(const prob_linear_doubling&) = delete;
        prob_linear_doubling(prob_linear_doubling&& rhs)  = default;
        prob_linear_doubling& operator=(prob_linear_doubling&& ) = default;

    private:
        using base_type::capacity;
        using base_type::thresh;
        using base_type::table;

        size_type bitmask;

        // Access Functions ********************************************************
        inline size_type index(size_type i) const { return i & bitmask; }
        inline size_type mod(size_type i)   const { return i & bitmask; }

        // Helper Function *********************************************************
        static size_type next_power_of_two(size_type i)
        {
            size_type t = 2048;
            while (t < i) t <<= 1;
            return t;
        }

        // Growing *****************************************************************
        explicit prob_linear_doubling(size_type cap, size_type lthresh, this_type* )
            : base_type(cap, 1.), bitmask(cap-1)
        {
            thresh = lthresh;
        }

        inline void grow()
        {
            auto nsize  = capacity << 1;
            auto ntable = this_type(nsize, nsize*0.6, this);

            for (size_type i = 0; i <= bitmask; ++i)
            {
                auto temp = table[i];
                if (temp.first)
                {
                    ntable.insert(temp);
                }
            }

            (*this) = std::move(ntable);
        }
    };




    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_linear_doubling<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_linear_doubling<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };





/* Using Classic Linear Probing to Fill a Table Densely ***********************/
    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_linear : public prob_traits<prob_linear<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_linear<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;
    public:
        using size_type   = typename base_type::size_type;
        using key_type    = typename prob_traits<this_type>::key_type;
        using mapped_type = typename prob_traits<this_type>::mapped_type;

        prob_linear(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
            : base_type(std::max<size_type>(cap, 500), size_constraint)
        {
            //acap = capacity-300;
        }
        prob_linear(const prob_linear&) = delete;
        prob_linear& operator=(const prob_linear&) = delete;
        prob_linear(prob_linear&& rhs)  = default;
        prob_linear& operator=(prob_linear&& ) = default;

    private:
        using base_type::alpha;
        using base_type::n;
        using base_type::capacity;
        using base_type::table;
        using base_type::h;

        // size_type acap;

        static constexpr size_type bitmask = (1ull<<32)-1;


        // Access Functions ********************************************************
        inline size_type index(size_type i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_type mod(size_type i)   const
        { return (i < capacity) ? i : i-capacity; }

        // Growing *****************************************************************
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

        // Specialized Deletion stuff **********************************************

        inline void propagate_remove(const size_type hole)
        {
            size_type thole = hole;
            for (size_type i = hole+1; ; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 ) break;
                auto tind = h(temp.first);
                if (tind <= thole)
                {
                    table[thole] = temp;
                    thole = i;
                }
            }
            table[thole] = std::make_pair(0,0);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_linear<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_linear<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };










// *****************************************************************************
// Same as Above, but Growing Using in Place Migration *************************
// *****************************************************************************

    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_linear_inplace : public prob_traits<prob_linear_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_linear_inplace<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;

    public:
        using size_type   = typename base_type::size_type;
        using key_type    = typename prob_traits<this_type>::key_type;
        using mapped_type = typename prob_traits<this_type>::mapped_type;
    private:
        using value_intern = std::pair<key_type, mapped_type>;

        static constexpr size_type max_size = 16ull << 30;

    public:
        prob_linear_inplace(size_type cap = 0, double size_constraint = 1.1, size_type /*steps*/=0)
            : base_type(0, size_constraint), bla(0)
        {
            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            if (temp) table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            thresh   = (cap) ? cap*beta  : 2048*beta;
            //acap = capacity-300;

            std::fill(table.get(), table.get()+capacity, value_intern());
        }
        prob_linear_inplace(const prob_linear_inplace&) = delete;
        prob_linear_inplace& operator=(const prob_linear_inplace&) = delete;
        prob_linear_inplace(prob_linear_inplace&& rhs)  = default;
        prob_linear_inplace& operator=(prob_linear_inplace&& ) = default;

    private:
        using base_type::alpha;
        using base_type::beta;
        using base_type::n;
        using base_type::capacity;
        using base_type::thresh;
        using base_type::table;
        using base_type::dec_n;

        //size_type acap;
        //double    factor;
        size_type bla;

        static constexpr size_type bitmask = (1ull<<32)-1;

        using base_type::h;

        // Access Functions ********************************************************
        inline size_type index(size_type i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_type mod(size_type i)   const
        { return (i < capacity) ? i : i - capacity; }

    public:
        using base_type::insert;

    private:
        // Growing *****************************************************************
        inline void grow()
        {
            // auto ntable = this_type(n, alpha);

            size_type ncap    = n*alpha;
            size_type nthresh = n*beta;
            //double    nfactor = double(ncap-300)/double(1ull << 32);

            std::fill(table.get()+capacity, table.get()+ncap, value_intern());

            auto old_cap = capacity;
            capacity = ncap;
            thresh   = nthresh;
            //acap     = ncap-300;
            n = 0;

            std::vector<value_intern> buffer;

            for (int i = old_cap - 1; i >= 0; --i)
            {
                auto temp = table[i];

                if (temp.first)
                {
                    table[i] = value_intern();
                    auto ind = h(temp.first);
                    if (ind >= size_t(i))
                        insert(temp);
                    else
                        buffer.push_back(temp);
                }
                else if (! buffer.empty())
                {
                    bla = std::max(bla, buffer.size());
                    for (auto it = buffer.begin(); it != buffer.end(); it++)
                    {
                        insert(*it);
                    }
                    buffer.clear();
                }
            }

            if (! buffer.empty())
            {
                bla = std::max(bla, buffer.size());
                for (auto it = buffer.begin(); it != buffer.end(); it++)
                {
                    insert(*it);
                }
            }
        }

        // Specialized Deletion stuff **********************************************

        inline void propagate_remove(const size_type hole)
        {
            size_type thole = hole;
            for (size_type i = hole+1; ; ++i)
            {
                if (i == capacity)
                {
                    table[thole] = std::make_pair(0,0);
                    // wrap around is hard and uncommon
                    // thus do something bad we reinsert elements from start to ...
                    for (int i = 0; ; ++i)
                    {
                        auto temp = table[i];
                        table[i] = std::make_pair(0,0);
                        dec_n();
                        insert(temp);
                    }
                    return;
                }
                auto temp = table[i];

                if ( temp.first == 0 ) break;
                auto tind = h(temp.first);
                if (tind <= thole)
                {
                    table[thole] = temp;
                    thole = i;
                }
            }
            table[thole] = std::make_pair(0,0);
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(8) << "busize";
            base_type::print_init_header(out);
        }

        inline        void print_init_data  (otm::output_type& out)
        {
            out << otm::width(8) << bla;
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_linear_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_linear_inplace<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };


} // namespace dysect
