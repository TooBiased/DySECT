#pragma once

/*******************************************************************************
 * include/prob_robin.h
 *
 * robin_prob and robin_prob_inplace implement robin hood hashing, a
 * variant of linear probing.  The inplace variant uses memory
 * overcommiting to resize the table without full table reallocations,
 * that would temporarily violate the memory constraint.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.hpp"
#include "utils/output.hpp"

#include "prob_base.h"

namespace otm = utils_tm::out_tm;

namespace dysect
{

    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = triv_config>
    class prob_robin : public prob_traits<prob_robin<K,D,HF,Conf> >::base_type
    {
    private:
        using  this_type = prob_robin<K,D,HF,Conf>;
        using  base_type = typename prob_traits<this_type>::base_type;

        friend base_type;
    public:
        using size_type      = typename base_type::size_type;
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

        prob_robin(size_type cap = 0      , double size_constraint = 1.1,
                  size_type /*dis_steps*/ = 0, size_type /*seed*/ = 0)
            : base_type(std::max<size_type>(cap, 500), size_constraint),
              pdistance(0)
        { factor = double(capacity-300)/double(1ull << 32); }

        prob_robin(const prob_robin&) = delete;
        prob_robin& operator=(const prob_robin&) = delete;

        prob_robin(prob_robin&& rhs)  = default;
        prob_robin& operator=(prob_robin&& ) = default;

    private:
        using base_type::alpha;
        using base_type::n;
        using base_type::capacity;
        using base_type::hasher;
        using base_type::table;

        size_type pdistance;
        double factor;

        static constexpr size_type bitmask = (1ull << 32) - 1;

        using base_type::h;
        using base_type::inc_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

    public:

        //specialized functions because of Robin Hood Hashing
        inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
        {
            return insert(std::make_pair(k,d));
        }

        inline std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>& t)
        {
            // using doubles makes the element order independent from the capacity
            // thus growing gets even easier
            double ind     = dindex(hasher(t.first));
            auto   current = t;

            for (size_type i = ind; ; ++i)
            {
                auto temp  = table[i];

                if ( temp.first == t.first )
                {
                    return std::make_pair(make_iterator(&table[i]), false);
                }
                if ( temp.first == 0 )
                {
                    if (i == capacity - 1)
                        return std::make_pair(base_type::end(), false);
                    table[i] = current;
                    inc_n();
                    pdistance = std::max<size_type>(pdistance, i-size_type(ind));
                    return std::make_pair(make_iterator(&table[i]), true);
                }
                double tind = dindex(hasher(temp.first));
                if ( tind > ind )
                {
                    std::swap(table[i], current);
                    pdistance = std::max<int>(pdistance, i-size_type(ind));
                    ind = tind;
                }
            }
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    return make_iterator(&table[i]);
                }
            }
            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    return make_citerator(&table[i]);
                }
            }
            return base_type::cend();
        }

        inline size_type erase(const key_type& k)
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance ; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    base_type::dec_n();
                    propagate_remove(i);
                    return 1;
                }
            }
            return 0;
        }

        inline int displacement(const key_type& k) const
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    return i-ind;
                }
            }
            return -1;
        }

    private:
        inline size_type index (size_type i) const
        { return double(bitmask & i) * factor; }
        inline double dindex(size_type i) const
        { return double(bitmask & i) * factor; }
        inline size_type mod(size_type i)    const
        { return i; }


        inline void grow()
        {
            auto ntable = this_type(n, alpha);

            size_type tn = n;
            size_type tdistance = ntable.migrate(*this);

            (*this) = std::move(ntable);

            n         = tn;
            pdistance = tdistance;
        }

        inline size_type migrate(this_type& source)
        {
            size_type distance   = 0;
            size_type target_pos = 0;
            for (size_type i = 0; i < source.capacity; ++i)
            {
                auto current = source.table[i];
                if (!current.first) continue;
                auto hash = h(current.first);
                if (target_pos > hash)
                    distance   = std::max<size_type>(distance, target_pos - hash);
                else
                    target_pos = hash;
                table[target_pos++] = current;
            }
            return distance;
        }

        inline void propagate_remove(const size_type hole)
        {
            size_type thole = hole;
            for (size_type i = hole+1; ; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 ) break;
                auto tind = h(temp.first);
                if ( tind >= i)       break;

                table[thole] = temp;
                thole = i;
            }
            table[thole] = std::make_pair(0,0);
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(6) << "pdis";
            base_type::print_init_header(out);
        }

        inline        void print_init_data  (otm::output_type& out)
        {
            out << otm::width(6) << pdistance;
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_robin<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_robin<K,D,HF,Conf>;
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
    class prob_robin_inplace : public prob_traits<prob_robin_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_robin_inplace<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        friend base_type;

    public:
        using size_type      = typename base_type::size_type;
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;
    private:
        using value_intern   = typename base_type::value_intern;
        static constexpr size_type max_size = 16ull << 30;
    public:

        prob_robin_inplace(size_type cap = 0      , double size_constraint = 1.1,
                       size_type /*dis_steps*/ = 0, size_type /*seed*/ = 0)
            : base_type(0, size_constraint),
              pdistance(0)
        {
            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            if (temp) table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            thresh   = (cap) ? cap*beta  : 2048*beta;
            factor = double(capacity-300)/double(1ull << 32);

            std::fill(table.get(), table.get()+capacity, value_intern());
        }

        prob_robin_inplace(const prob_robin_inplace&) = delete;
        prob_robin_inplace& operator=(const prob_robin_inplace&) = delete;

        prob_robin_inplace(prob_robin_inplace&& rhs)  = default;
        prob_robin_inplace& operator=(prob_robin_inplace&& ) = default;

    private:
        using base_type::alpha;
        using base_type::beta;
        using base_type::n;
        using base_type::capacity;
        using base_type::thresh;
        using base_type::hasher;
        using base_type::table;

        size_type pdistance;
        double factor;

        static constexpr size_type bitmask = (1ull << 32) - 1;

        using base_type::h;
        using base_type::inc_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

    public:

        //specialized functions because of Robin Hood Hashing
        inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
        {
            return insert(std::make_pair(k,d));
        }

        inline std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>& t)
        {
            // using doubles makes the element order independent from the capacity
            // thus growing gets even easier
            double ind     = dindex(hasher(t.first));
            auto   current = t;

            for (size_type i = ind; ; ++i)
            {
                auto temp  = table[i];

                if ( temp.first == t.first )
                {
                    return std::make_pair(make_iterator(&table[i]), false);
                }
                if ( temp.first == 0 )
                {
                    if (i == capacity - 1)
                        return std::make_pair(base_type::end(), false);
                    table[i] = current;
                    inc_n();
                    pdistance = std::max<size_type>(pdistance, i-size_type(ind));
                    return std::make_pair(make_iterator(&table[i]), true);
                }
                double tind = dindex(hasher(temp.first));
                if ( tind > ind )
                {
                    std::swap(table[i], current);
                    pdistance = std::max<int>(pdistance, i-size_type(ind));
                    ind = tind;
                }
            }
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    return make_iterator(&table[i]);
                }
            }
            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    return make_citerator(&table[i]);
                }
            }
            return base_type::cend();
        }

        inline size_type erase(const key_type& k)
        {
            auto ind = h(k);

            for (size_type i = ind; i <= ind+pdistance ; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 )
                {
                    break;
                }
                else if ( temp.first == k )
                {
                    base_type::dec_n();
                    propagate_remove(i);
                    return 1;
                }
            }
            return 0;
        }

    private:
        inline size_type index (size_type i) const
        { return double(bitmask & i) * factor; }
        inline double dindex(size_type i) const
        { return double(bitmask & i) * factor; }
        inline size_type mod(size_type i)    const
        { return i; }

        inline void grow()
        {
            size_type ncap    = n*alpha;
            size_type nthresh = n*beta;
            double    nfactor = double(ncap-300)/double(1ull << 32);

            std::fill(table.get()+capacity, table.get()+ncap, value_intern());

            size_type ocap = capacity;
            capacity = ncap;
            thresh   = nthresh;
            factor   = nfactor;

            size_type tn = n;
            n = 0;

            migrate(ocap);

            n = tn;
        }

        inline void migrate(size_type index)
        {
            std::vector<value_intern> buffer;

            size_type distance   = 0;
            for (int i = index; i >= 0; --i)
            {
                auto temp = table[i];

                if (!temp.first)
                    continue;

                table[i] = value_intern();

                auto nind = h(temp.first);
                if (nind < size_type(i)) { buffer.push_back(temp); continue; }

                size_type t = nind;

                while (temp.first)
                {
                    std::swap(table[t++], temp);
                }
                pdistance = std::max<size_type>(t-i, distance);
            }
            for (auto it = buffer.begin(); it != buffer.end(); it++)
            {
                insert(*it);
            }

        }

        inline void propagate_remove(const size_type hole)
        {
            size_type thole = hole;
            for (size_type i = hole+1; ; ++i)
            {
                auto temp = table[i];

                if ( temp.first == 0 ) break;
                auto tind = h(temp.first);
                if ( tind >= i)       break;

                table[thole] = temp;
                thole = i;
            }
            table[thole] = std::make_pair(0,0);
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(6) << "pdis";
            base_type::print_init_header(out);
        }

        inline        void print_init_data  (otm::output_type& out)
        {
            out << otm::width(6) << pdistance;
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_robin_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_robin_inplace<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };

} // namespace dysect
