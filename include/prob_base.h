#pragma once

/*******************************************************************************
 * include/prob_base.h
 *
 * prob_base is similar to cuckoo_base, in that it encapsules
 * everything, that all probing based hash tables have in common.
 *
 * Part of Project DySECT - https://github.com/TooBiased/DySECT.git
 *
 * Copyright (C) 2017 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <functional>
#include <memory>
#include <vector>
#include <tuple>

#include "utils/fastrange.hpp"
#include "utils/output.hpp"

#include "bucket.h"
#include "iterator_base.h"

namespace otm = utils_tm::out_tm;

namespace dysect
{
    struct triv_config
    { };

    template <class T>
    class prob_traits;
    template <class T>
    class iterator_incr;





    template<class SpProb>
    class prob_base
    {
    private:
        using this_type         = prob_base<SpProb>;
        using specialized_type  = typename prob_traits<SpProb>::specialized_type;
        using hash_function_type      = typename prob_traits<SpProb>::hash_function_type;

        friend specialized_type;
        friend iterator_incr<this_type>;

    public:
        using size_type      = size_t;
        using key_type       = typename prob_traits<SpProb>::key_type;
        using mapped_type    = typename prob_traits<SpProb>::mapped_type;
        using iterator       = iterator_base<iterator_incr<this_type> >;
        using const_iterator = iterator_base<iterator_incr<this_type>, true>;

    private:
        using value_intern         = std::pair<key_type,mapped_type>;

    public:
        prob_base(size_type cap, double alpha)
            : alpha(alpha), beta((alpha+1.)/2.), n(0),
              capacity((cap) ? cap*alpha : 2048*alpha),
              thresh  ((cap) ? cap*beta  : 2048*beta)
        {
            if (cap) table = std::make_unique<value_intern[]>(capacity);
        }

        ~prob_base() = default;

        prob_base(const prob_base&) = delete;
        prob_base& operator=(const prob_base&) = delete;

        prob_base(prob_base&& rhs) = default;
        prob_base& operator=(prob_base&& rhs) = default;

    private:
        /*** members that should become private at some point *********************/
        double     alpha;
        double     beta;
        size_type  n;
        size_type  capacity;
        size_type  thresh;
        hash_function_type  hasher;

        std::unique_ptr<value_intern[]> table;

    public:
        // Basic Hash Table Functionality ******************************************
        iterator                  find  (const key_type& k);
        const_iterator            find  (const key_type& k) const;
        std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
        std::pair<iterator, bool> insert(const value_intern& t);
        size_type                 erase (const key_type& k);
        int                 displacement(const key_type& k) const;

        // Easy use Accessors for std compliance ***********************************
        inline iterator           begin ()
        {
            auto temp = make_iterator(&table[0]);
            if (!temp->first) temp++;
            return temp;
        }
        inline const_iterator     begin()  const
        {
            return static_cast<specialized_type*>(this)->cbegin();
        }
        inline const_iterator     cbegin() const
        {
            auto temp = make_citerator(&table[0]);
            if (!temp->first) temp++;
            return temp;
        }
        inline iterator           end   () { return make_iterator(nullptr);  }
        inline const_iterator     end   () const { return static_cast<specialized_type*>(this)->cend(); }
        inline const_iterator     cend  () const { return make_citerator(nullptr); }

        mapped_type&              at    (const key_type& k);
        const mapped_type&        at    (const key_type& k) const;
        mapped_type&              operator[](const key_type& k);
        size_type                 count (const key_type& k) const;

        size_type                 get_capacity() const { return capacity; }

    private:
        // Easy iterators **********************************************************
        inline iterator make_iterator(value_intern* pair) const
        { return iterator(pair, *this); }
        inline const_iterator make_citerator(value_intern* pair) const
        { return const_iterator(pair, *this); }

        // implementation specific functions (static polymorph) ********************
        inline size_type h(key_type k) const
        { return static_cast<const specialized_type*>(this)->index(hasher(k)); }

        inline void inc_n()
        { if (++n > thresh) static_cast<specialized_type*>(this)->grow(); }
        inline void dec_n() { --n; }

        // Private helper function *************************************************
        void propagate_remove(size_type origin);

    public:
        inline static void print_init_header(otm::output_type& out)
        { out << otm::width(10) << "f_cap"; }
        inline void print_init_data  (otm::output_type& out)
        { out << otm::width(10) << capacity;}
    };



// Implementation of main functionality ****************************************

    template<class SpProb>
    inline typename prob_base<SpProb>::iterator
    prob_base<SpProb>::find(const key_type& k)
    {
        auto ind = h(k);

        for (size_type i = ind; ; ++i)
        {
            size_type ti = static_cast<specialized_type*>(this)->mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 )
            {
                break;
            }
            else if ( temp.first == k )
            {
                return make_iterator(&table[ti]);
            }
        }
        return end();
    }

    template<class SpProb>
    inline typename prob_base<SpProb>::const_iterator
    prob_base<SpProb>::find(const key_type& k) const
    {
        auto ind = h(k);

        for (size_type i = ind; ; ++i)
        {
            size_type ti = static_cast<const specialized_type*>(this)->mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 )
            {
                break;
            }
            else if ( temp.first == k )
            {
                return make_citerator(&table[ti]);
            }
        }
        return cend();
    }

    template<class SpProb>
    inline std::pair<typename prob_base<SpProb>::iterator, bool>
    prob_base<SpProb>::insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    template<class SpProb>
    inline std::pair<typename prob_base<SpProb>::iterator, bool>
    prob_base<SpProb>::insert(const value_intern& t)
    {
        auto ind = h(t.first);

        for (size_type i = ind; ; ++i)
        {
            size_type ti = static_cast<specialized_type*>(this)->mod(i);
            auto temp = table[ti];
            if ( temp.first == t.first)
            {
                return std::make_pair(make_iterator(&table[ti]), false);
            }
            if ( temp.first == 0 )
            {
                table[ti] = t;
                //hcounter.add(i - ind);
                static_cast<SpProb*>(this)->inc_n();
                return std::make_pair(make_iterator(&table[ti]), true);
            }
        }
        return std::make_pair(end(), false);
    }

    template<class SpProb>
    inline typename prob_base<SpProb>::size_type
    prob_base<SpProb>::erase(const key_type& k)
    {
        auto ind = h(k);

        for (size_type i = ind; ; ++i)
        {
            size_type ti = static_cast<specialized_type*>(this)->mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 )
            {
                break;
            }
            else if ( temp.first == k )
            {
                dec_n();
                table[ti] = std::make_pair(0,0);
                static_cast<SpProb*>(this)->propagate_remove(ti);
                return 1;
            }
        }
        return 0;
    }

    template<class SpProb>
    inline int
    prob_base<SpProb>::displacement(const key_type& k) const
    {
        auto ind = h(k);

        for (size_type i = ind; ; ++i)
        {
            size_type ti = static_cast<const specialized_type*>(this)->mod(i);
            auto temp = table[ti];

            if ( temp.first == 0 )
            {
                break;
            }
            else if ( temp.first == k )
            {
                return i - ind;
            }
        }
        return -1;
    }


// Accessor Implementations ****************************************************

    template<class SpProb>
    inline typename prob_base<SpProb>::mapped_type&
    prob_base<SpProb>::at(const key_type& k)
    {
        auto a = static_cast<specialized_type*>(this)->find(k);
        if (a == end()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class SpProb>
    inline const typename prob_base<SpProb>::mapped_type&
    prob_base<SpProb>::at(const key_type& k) const
    {
        auto a = static_cast<const specialized_type*>(this)->find(k);
        if (a == cend()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class SpProb>
    inline typename prob_base<SpProb>::mapped_type&
    prob_base<SpProb>::operator[](const key_type& k)
    {
        auto t = static_cast<specialized_type*>(this)->insert(k, mapped_type());
        return (*t.first).second;
    }

    template<class SpProb>
    inline typename prob_base<SpProb>::size_type
    prob_base<SpProb>::count(const key_type& k) const
    {
        return (static_cast<const specialized_type*>(this)->find(k) != cend()) ? 1 : 0;
    }



// Private Help Function *******************************************************

    template<class SpProb>
    inline void prob_base<SpProb>::propagate_remove(const size_type origin)
    {
        size_type tempn = n;
        n = 0;
        for (size_type i = origin+1; ; ++i)
        {
            size_type ti = static_cast<const SpProb*>(this)->mod(i);
            auto temp = table[ti];

            if (temp.first == 0) break;

            table[ti] = std::make_pair(0,0);
            insert(temp);
        }
        n = tempn;
    }



// Iterator increment **********************************************************

    template<class Specialized>
    class iterator_incr<prob_base<Specialized> >
    {
    public:
        using table_type  = prob_base<Specialized>;

    private:
        using key_type    = typename table_type::key_type;
        using mapped_type = typename table_type::mapped_type;
        using ipointer    = std::pair<const key_type, mapped_type>*;

    public:
        iterator_incr(const table_type& table_)
            : end_ptr(reinterpret_cast<ipointer>
                      (&table_.table[table_.capacity - 1]))
        { }
        iterator_incr(const iterator_incr&) = default;
        iterator_incr& operator=(const iterator_incr&) = default;

        ipointer next(ipointer cur)
        {
            while (cur < end_ptr)
            {
                if ((++cur)->first) return cur;
            }
            return nullptr;
        }

    private:
        const ipointer end_ptr;
    };

}
