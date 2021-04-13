#pragma once

/*******************************************************************************
 * include/prob_hops.h
 *
 * hops_prob and hops_prob_inplace implement the hopscotch hashing
 * scheme.  The inplace variant uses memory overcommiting to resize
 * the table without full table reallocations, that would temporarily
 * violate the memory constraint.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "utils/default_hash.h"
#include "prob_base.h"

namespace dysect
{

    template<class IntegerType = uint>
    struct coalesced_config
    {
        using integer_type = IntegerType;
    };







// *****************************************************************************
// MAIN CLASS ******************************************************************
// *****************************************************************************

    template <class K, class D, class HF = hash::default_hash,
              class Conf = coalesced_config<> >
    class prob_coalesced : public prob_traits<prob_coalesced<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type   = prob_coalesced<K,D,HF,Conf>;
        using base_type   = typename prob_traits<this_type>::base_type;
        using offset_type = typename Conf::integer_type;

        friend base_type;

    public:
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

        prob_coalesced(size_t cap = 0      , double size_constraint = 1.1,
                 size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
            : base_type(cap, size_constraint)
        {
            offset_table = std::make_unique<integer_type[]>(capacity);
            std::fill(offset_table.get(), offset_table.get() + capacity, 0);
        }

        prob_coalesced(const prob_coalesced&) = delete;
        prob_coalesced& operator=(const prob_coalesced&) = delete;

        prob_coalesced(prob_coalesced&& rhs)  = default;
        prob_coalesced& operator=(prob_coalesced&& ) = default;

    private:
        using base_type::alpha;
        using base_type::n;
        using base_type::capacity;
        using base_type::table;

        std::unique_ptr<offset_type[]> offset_table;

        static constexpr size_t bitmask = (1ull << 32) - 1;

        using base_type::h;
        using base_type::inc_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        //specialized functions because of Hops Hashing
        inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
        {
            return insert(std::make_pair(k,d));
        }

        inline std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>& t)
        {
            auto ind = h(k);
            auto current = table[ind];
            if (ind != h(current.first))
            {
                // current slot is used by another queue
                size_t eind = ind;
                for (; ; eind = mod(eind+1))
                {
                    if (!table[eind].first) break;
                }
                auto pind = h(table[ind].first);
                auto next = pind;
                do
                {
                    pind = next;
                    next = pind + offset_table[pind];
                } while (next != ind);
                offset_table[pind] = eind - pind;

                table[eind]        = current;
                offset_table[eind] = offset_table[ind];
                table[ind]         = t;
                offset_table[ind]  = 0;
            }

            do
            {
                if (temp.first == k)
                    return std::make_pair(make_iterator(&table[ind]), false);

                auto offs = offset_table[ind];
                ind += offs;
                current = table[ind]
            } while (offs);

            for (size_t i=mod(ind+1); ;i+=mod(i+1))
            {
                auto temp = table[i];

                if (temp.first == 0)
                {
                    table[i] = t;
                    // if (i-ind > std::numeric_limits<offset_type>::max())
                    //     move some elements around
                    offset_table[ind] = i-ind;
                    return std::make_pair(make_iterator(&table[ind]), true);
                }
            }

            // should be unreachable code
            // (a full table leads to an infinite loop)
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            auto ind = h(k);
            auto current = table[ind];
            if (ind != h(current.first)) return end();

            do
            {
                if (temp.first == k)
                    return make_iterator(&table[ind]);

                auto offs = offset_table[ind];
                ind += offs;
                current = table[ind];
            } while (offs);

            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            auto ind = h(k);
            auto current = table[ind];
            if (ind != h(current.first)) return cend();

            do
            {
                if (temp.first == k)
                    return make_citerator(&table[ind]);

                auto offs = offset_table[ind];
                ind += offs;
                current = table[ind];
            } while (offs);

            return base_type::cend();
        }

        inline size_t erase(const key_type& k)
        {
            auto ind = h(k);
            auto current = table[ind];
            if (ind != h(current.first)) return cend();

            auto prev = ind;
            do
            {
                auto offs = offset_table[ind];
                if (temp.first == k)
                {
                    auto next = ind + offs;

                    if (prev == ind)
                    {
                        // we removed the first chain link
                        // switch elements with next and remove next
                        // still works for chain length 1 (prev == ind == next)
                        table[ind] = table[next];
                        ind = next;
                        offs = offset_table[next];
                        next = ind + offs;
                    }
                    offset_table[prev] = (offs) ? next - prev : 0;
                    table[ind] = std::make_pair(0,0);

                }

                prev = ind;
                mod(ind += offs);
                current = table[ind];
            } while (offs);

            return base_type::cend();


            auto ind = h(k);
            auto prev = ind;
            do
            {
                auto temp = table[ind];
                auto offs = offset_table[ind];

                if (temp.first == k)
                {
                    table[ind] = std::make_paier(0,0);
                    if (prev == ind)
                    {    // we deleted the first chain element
                        if (offs)
                        {
                            auto nexti = ind + offs;
                            auto nexte = table[nexti];
                            auto nexto = offset_table[nexti];
                            table[ind] = nexte;
                            table[off] = nexto+offs;
                        }
                    }
                    else if (offs == 0)
                    {
                        // we deleted the last chain element
                        offset_table[prev] = 0;
                    }
                    else
                    {
                        offset_table[ind] = 0;
                        auto next = ind + offs - prev;
                        offset_table[prev] = next;
                    }
                    return 1;
                }
                prev = ind;
                mod(ind += offs);

            } while (offs);

            return 0;
        }

    private:
        inline size_t index(size_t i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_t mod  (size_t i) const
        { return (i < capacity) ? i : i-capacity; }

        inline void grow()
        {
            auto ntable = this_type(n, alpha);

            for (size_t i = 0; i < capacity; ++i)
            {
                auto current  = table[i];
                if (current.first)
                    ntable.insert(current);
            }

            (*this) = std::move(ntable);
        }

    public:
        inline static void print_init_header(std::ostream& out)
        {
            base_type::print_init_header(out);
        }

        inline void print_init_data(std::ostream& out)
        {
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_coalesced<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_coalesced<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };











// *****************************************************************************
// Same as Above, but Growing Using in Place Migration *************************
// *****************************************************************************

    template <class K, class D, class HF = hash::default_hash,
              class Conf = coalesced_config<> >
    class prob_coalesced_inplace : public prob_traits<prob_coalesced_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_coalesced_inplace<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        static constexpr size_t nh_size = Conf::neighborhood_size;
        static constexpr size_t bitset_size = (nh_size < 64) ? nh_size : 64;
        static constexpr size_t bucket_size = (nh_size < 64) ? 1 : nh_size/64;
        using AugData_t = augment_data<bitset_size>;

        friend base_type;

    public:
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;
    private:
        using value_intern   = typename base_type::value_intern;

        static constexpr size_t max_size = 1ull << 30;

    public:
        prob_coalesced_inplace(size_t cap = 0      , double size_constraint = 1.1,
                        size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
            : base_type(0, size_constraint), nh_data(max_size, 0)
        {
            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            capacity = capacity/bucket_size;
            capacity = capacity*bucket_size;
            thresh   = (cap) ? cap*beta  : 2048*beta;

            std::fill(table.get(), table.get()+capacity, value_intern());
            nh_data.clear_init(capacity);
        }

        prob_coalesced_inplace(const prob_coalesced_inplace&) = delete;
        prob_coalesced_inplace& operator=(const prob_coalesced_inplace&) = delete;
        prob_coalesced_inplace(prob_coalesced_inplace&& rhs)  = default;
        prob_coalesced_inplace& operator=(prob_coalesced_inplace&& ) = default;

    private:
        using base_type::alpha;
        using base_type::beta;
        using base_type::n;
        using base_type::capacity;
        using base_type::thresh;
        using base_type::table;

        /* TODO implement offset table */

        static constexpr size_t bitmask = (1ull << 32) - 1;

        using base_type::h;
        using base_type::inc_n;
        using base_type::make_iterator;
        using base_type::make_citerator;

    public:
        //specialized functions because of Hops Hashing
        inline std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d)
        {
            return insert(std::make_pair(k,d));
        }

        inline std::pair<iterator, bool> insert(const std::pair<key_type, mapped_type>& t)
        {
            // we first have to check if t.first is already present
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            return base_type::cend();
        }

        inline size_t erase(const key_type& k)
        {
            return 0;
        }

    private:
        inline size_t index(size_t i) const
        { return utils_tm::fastrange64(capacity, i); }
        inline size_t mod  (size_t i) const
        { return (i < capacity) ? i : i-capacity; }

        inline void grow()
        {

            // TODO THIS METHOD IS UNTESTED
            // There are some weird corner cases, it looks alright though
            size_t osize = capacity;

            capacity = n*alpha;
            thresh   = n*beta;

            std::fill(table.get()+osize, table.get()+capacity, value_intern());
            // reset all offsets
            std::fill(offset_table.get(), offset_table.get()+capacity, 0);

            std::vector<value_intern> buffer;

            n = 0;

            for (int i = osize; i >= 0; --i)
            {
                auto current  = table[i];
                if (current.first)
                {
                    table[i] = value_intern();
                    if (h(current.first) > size_t(i))
                        insert(current);
                    else
                        buffer.push_back(current);
                }
            }
            for (auto it = buffer.begin(); it != buffer.end(); ++it)
            {
                insert(*it);
            }
        }

    public:
        inline static void print_init_header(std::ostream& out)
        {
            base_type::print_init_header(out);
        }

        inline void print_init_data(std::ostream& out)
        {
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_coalesced_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_coalesced_inplace<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };

} // namespace dysect
