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

#include "utils/default_hash.hpp"
#include "utils/output.hpp"

#include "prob_base.h"

namespace otm = utils_tm::out_tm;

namespace dysect
{

    template<size_t NS = 64>
    struct hopscotch_config
    {
        static constexpr size_t neighborhood_size = NS;
    };

    template <class AugmentData>
    class augment_data_accessor
    {
    public:
        static constexpr size_t bitmask = AugmentData::bitmask;
        static constexpr size_t nh_size = AugmentData::nh_size;

        augment_data_accessor(size_t* data) : data(*data) { }

        size_t& data;

        inline size_t get_neighborhood() const
        {
            return data & bitmask;
        }
        inline void set(size_t i)
        {
            data |= (1ull << i);
        }
        inline void unset(size_t i)
        {
            data &= ~(1ull << i);
        }
    };

    template <size_t ns>
    class augment_data
    {
        using uchar = uint8_t;
    public:
        static_assert(ns >   0, "Augment Data cannot handle neighborhood size 0!");
        static_assert(ns <= 64, "Augment Data cannot handle neighborhoods >64!");

        static constexpr size_t bitmask = (ns == 64) ? ~size_t(0) : (1ull << ns) -1;
        static constexpr size_t n_bytes = (ns>>3)+((ns&7) ? 1:0);
        static constexpr size_t nh_size = ns;
        using this_type     = augment_data<ns>;
        using augment_accessor_type = augment_data_accessor<this_type>;

        augment_data(size_t capacity, size_t initialized)
        {
            auto glob_n_bytes = capacity   *n_bytes+8-n_bytes;
            auto init_n_bytes = initialized*n_bytes+8-n_bytes;

            auto temp = static_cast<uchar*>(operator new (glob_n_bytes));
            data = std::unique_ptr<uchar[]>(temp);

            std::fill(data.get(), data.get()+init_n_bytes, 0);
            init = init_n_bytes;
        }

        augment_data(const augment_data&) = delete;
        augment_data& operator=(const augment_data&) = delete;

        augment_data(augment_data&& rhs) noexcept : data(nullptr)
        {   std::swap(data, rhs.data); }
        augment_data& operator=(augment_data&& rhs) noexcept
        {   std::swap(data, rhs.data); return *this; }

        inline augment_accessor_type get_accessor(size_t index)
        { return augment_accessor_type(reinterpret_cast<size_t*>(data.get() + index*n_bytes)); }

        inline const augment_accessor_type get_caccessor(size_t index) const
        { return augment_accessor_type(reinterpret_cast<size_t*>(data.get() + index*n_bytes)); }

        inline size_t get_neighborhood(size_t index) const
        {
            size_t temp = *reinterpret_cast<const size_t*>(data.get() + index*n_bytes);
            return temp & bitmask;
        }

        inline void clear_init(size_t upper)
        {
            auto up_n_bytes = upper*n_bytes+8-n_bytes;
            std::fill(data.get(), data.get()+up_n_bytes, 0);
            init = up_n_bytes;
        }

        size_t init;
        std::unique_ptr<uchar[]> data;
    };






// *****************************************************************************
// MAIN CLASS ******************************************************************
// *****************************************************************************

    template <class K, class D, class HF = utils_tm::hash_tm::default_hash,
              class Conf = hopscotch_config<> >
    class prob_hopscotch : public prob_traits<prob_hopscotch<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_hopscotch<K,D,HF,Conf>;
        using base_type = typename prob_traits<this_type>::base_type;

        static constexpr size_t nh_size = (Conf::neighborhood_size <= 64) ?
            Conf::neighborhood_size : 64;
        using AugData_t = augment_data<nh_size>;

        friend base_type;

    public:
        using key_type       = typename prob_traits<this_type>::key_type;
        using mapped_type    = typename prob_traits<this_type>::mapped_type;
        using iterator       = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

        prob_hopscotch(size_t cap = 0      , double size_constraint = 1.1,
                 size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
            : base_type(cap, size_constraint), nh_data(capacity-nh_size+1, capacity-nh_size+1)
        {
            //factor = double(capacity-nh_size)/double(1ull << 32);
            acap = capacity-nh_size;
        }

        prob_hopscotch(const prob_hopscotch&) = delete;
        prob_hopscotch& operator=(const prob_hopscotch&) = delete;

        prob_hopscotch(prob_hopscotch&& rhs)  = default;
        prob_hopscotch& operator=(prob_hopscotch&& ) = default;

    private:
        using base_type::alpha;
        using base_type::n;
        using base_type::capacity;
        using base_type::table;

        //double    factor;
        size_t acap;
        AugData_t nh_data;

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
            size_t ind  = h(t.first);
            auto   aug  = nh_data.get_accessor(ind);
            size_t bits = aug.get_neighborhood();

            for (size_t i = ind; bits; ++i, bits>>=1)
            {
                if (!(bits & 1)) continue;
                else
                {
                    auto temp = table[i];
                    if ( temp.first == t.first )
                    {
                        return std::make_pair(make_iterator(&table[i]), false);
                    }
                }
            }

            for (size_t i = ind; ; ++i)
            {
                auto temp  = table[i];
                if ( temp.first == 0)
                {
                    size_t ti = i;
                    if (ti >= ind + nh_size)
                    {
                        bool successful;
                        std::tie(successful, ti) = move_gap(i, ind+nh_size);
                        if (!successful) break;
                    }
                    table[ti] = t;
                    aug.set(ti-ind);
                    inc_n();
                    return std::make_pair(make_iterator(&table[ti]), true);
                }
            }
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; ++i, bits>>=1)
            {
                if (!(bits & 1)) continue;
                auto temp = table[i];
                if ( temp.first == k )
                {
                    return make_iterator(&table[i]);
                }
            }
            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; ++i, bits>>=1)
            {
                if (!(bits & 1)) continue;
                auto temp = table[i];
                if ( temp.first == k )
                {
                    return make_citerator(&table[i]);
                }
            }
            return base_type::cend();
        }

        inline size_t erase(const key_type& k)
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; ++i, bits>>=1)
            {
                if (!(bits&1)) continue;
                auto tempk = table[i].first;
                if ( tempk == k )
                {
                    nh_data.get_accessor(ind).unset(i-ind);
                    table[i] = std::make_pair(0,0);
                    return 1;
                }
            }
            return 0;
        }

        inline int displacement(const key_type& k) const
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; ++i, bits>>=1)
            {
                if (!(bits & 1)) continue;
                auto temp = table[i];
                if ( temp.first == k )
                {
                    return i-ind;
                }
            }
            return -1;
        }

    private:
        inline size_t index(size_t i) const
        { return utils_tm::fastrange64(acap, i); }
        inline size_t mod  (size_t i) const
        { return i; }

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

        inline std::pair<bool, size_t> move_gap(const size_t pos, const size_t goal)
        {
            for (size_t i = pos - nh_size + 1; i < pos; ++i)
            {
                auto current = table[i];
                auto ind = h(current.first);
                if (ind + nh_size > pos)
                {
                    auto aug = nh_data.get_accessor(ind);
                    aug.unset(i  -ind);
                    aug.set  (pos-ind);

                    table[pos] = current;
                    table[i]   = std::make_pair(0,0);
                    if (i < goal) return std::make_pair(true, i);
                    else return move_gap(i, goal);
                }
            }
            return std::make_pair(false, pos);
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(6) << "nghb";
            base_type::print_init_header(out);
        }

        inline void print_init_data(otm::output_type& out)
        {
            out << otm::width(6) << nh_size;
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_hopscotch<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_hopscotch<K,D,HF,Conf>;
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
              class Conf = hopscotch_config<> >
    class prob_hopscotch_inplace : public prob_traits<prob_hopscotch_inplace<K,D,HF,Conf> >::base_type
    {
    private:
        using this_type = prob_hopscotch_inplace<K,D,HF,Conf>;
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
        prob_hopscotch_inplace(size_t cap = 0      , double size_constraint = 1.1,
                        size_t /*dis_steps*/ = 0, size_t /*seed*/ = 0)
            : base_type(0, size_constraint), nh_data(max_size, 0)
        {
            value_intern* temp = reinterpret_cast<value_intern*>(operator new (max_size));
            table = std::unique_ptr<value_intern[]>(temp);

            capacity = (cap) ? cap*alpha : 2048*alpha;
            capacity = capacity/bucket_size;
            capacity = capacity*bucket_size;
            thresh   = (cap) ? cap*beta  : 2048*beta;

            //factor = double(capacity/bucket_size-bitset_size)/double(1ull << 32);
            acap = capacity/bucket_size-bitset_size;

            std::fill(table.get(), table.get()+capacity, value_intern());
            nh_data.clear_init(capacity);
        }

        prob_hopscotch_inplace(const prob_hopscotch_inplace&) = delete;
        prob_hopscotch_inplace& operator=(const prob_hopscotch_inplace&) = delete;
        prob_hopscotch_inplace(prob_hopscotch_inplace&& rhs)  = default;
        prob_hopscotch_inplace& operator=(prob_hopscotch_inplace&& ) = default;

    private:
        using base_type::alpha;
        using base_type::beta;
        using base_type::n;
        using base_type::capacity;
        using base_type::thresh;
        using base_type::table;

        //double    factor;
        size_t    acap;
        AugData_t nh_data;

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
            size_t ind  = h(t.first);
            auto   aug  = nh_data.get_accessor(ind);
            size_t bits = aug.get_neighborhood();

            for (size_t i = ind; bits; i+=bucket_size, bits>>=1)
            {
                if (!(bits & 1)) continue;
                else
                {
                    for (size_t ti = 0; ti<bucket_size; ++ti)
                    {
                        auto temp = table[i+ti];
                        if ( temp.first == t.first )
                        {
                            return std::make_pair(make_iterator(&table[i]), false);
                        }
                    }
                }
            }

            for (size_t i = ind; ; ++i)
            {
                auto temp  = table[i];
                if ( temp.first == 0)
                {
                    size_t ti = i;
                    if (ti >= ind + nh_size)
                    {
                        bool successful;
                        std::tie(successful, ti) = move_gap(i, ind+nh_size);
                        if (!successful) break;
                    }
                    table[ti] = t;
                    aug.set((ti-ind)/bucket_size);
                    inc_n();
                    return std::make_pair(make_iterator(&table[ti]), true);
                }
            }
            return std::make_pair(base_type::end(), false);
        }

        inline iterator find(const key_type& k)
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; i+=bucket_size, bits>>=1)
            {
                if (!(bits & 1)) continue;
                for (size_t ti = 0; ti < bucket_size; ++ti)
                {
                    auto temp = table[i+ti];
                    if ( temp.first == k )
                    {
                        return make_iterator(&table[i+ti]);
                    }
                }
            }
            return base_type::end();
        }

        inline const_iterator find(const key_type& k) const
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; i+=bucket_size, bits>>=1)
            {
                if (!(bits & 1)) continue;
                for (size_t ti = 0; ti < bucket_size; ++ti)
                {
                    auto temp = table[i+ti];
                    if ( temp.first == k )
                    {
                        return make_citerator(&table[i+ti]);
                    }
                }
            }
            return base_type::cend();
        }

        inline size_t erase(const key_type& k)
        {
            auto ind = h(k);
            size_t bits = nh_data.get_neighborhood(ind);

            for (size_t i = ind; bits; i+=bucket_size, bits>>=1)
            {
                if (!(bits&1)) continue;
                for (size_t ti = 0; ti < bucket_size; ++ti)
                {
                    auto tempk = table[i + ti].first;
                    if ( tempk == k )
                    {
                        table[i] = std::make_pair(0,0);
                        for (size_t tti = 0; tti < bucket_size; ++tti)
                        {
                            auto ttk = table[i+tti].first;
                            if (h(ttk) == ind) return 1;
                        }
                        nh_data.get_accessor(ind).unset((i-ind)/bucket_size);
                        return 1;
                    }
                }
            }
            return 0;
        }

    private:
        inline size_t index(size_t i) const
        { return utils_tm::fastrange64(acap,i) * bucket_size; }
        inline size_t mod  (size_t i) const
        { return i; }

        inline void grow()
        {
            size_t osize = capacity;

            capacity = n*alpha;
            capacity = capacity/bucket_size;
            capacity = capacity*bucket_size;
            thresh   = n*beta;

            //factor = double(capacity/bucket_size-bitset_size)/double(1ull << 32);
            acap = capacity/bucket_size - bitset_size;

            std::fill(table.get()+osize, table.get()+capacity, value_intern());
            nh_data.clear_init(capacity);

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

        inline std::pair<bool, size_t> move_gap(const size_t pos, const size_t goal)
        {
            for (size_t i = pos - nh_size + 1; i < pos; ++i)
            {
                auto current = table[i];
                auto ind = h(current.first);
                if (ind + nh_size > pos)
                {
                    auto aug = nh_data.get_accessor(ind);
                    table[i] = std::make_pair(0,0);

                    // Make sure, the same bucket does not store another element
                    // hashed to ind! Then delete the bit
                    // This is not technically correct, since it might look in the
                    // next bucket, but faster this way, and additional bits don't
                    // hurt.
                    bool other = false;
                    for (size_t ti = 0; ti < bucket_size; ++ti)
                    {
                        if (h(table[i+ti].first) == ind) other = true;
                    }
                    if (! other)
                        aug.unset((i  -ind)/bucket_size);

                    aug.set  ((pos-ind)/bucket_size);

                    table[pos] = current;
                    if (i < goal) return std::make_pair(true, i);
                    else return move_gap(i, goal);
                }
            }
            return std::make_pair(false, pos);
        }

    public:
        inline static void print_init_header(otm::output_type& out)
        {
            out << otm::width(6) << "nghb";
            base_type::print_init_header(out);
        }

        inline void print_init_data(otm::output_type& out)
        {
            out << otm::width(6) << nh_size;
            base_type::print_init_data(out);
        }
    };


    template<class K, class D, class HF, class Conf>
    class prob_traits<prob_hopscotch_inplace<K,D,HF,Conf> >
    {
    public:
        using specialized_type   = prob_hopscotch_inplace<K,D,HF,Conf>;
        using base_type          = prob_base<specialized_type>;
        using hash_function_type = HF;
        using config_type        = Conf;

        using key_type      = K;
        using mapped_type   = D;
    };

} // namespace dysect
