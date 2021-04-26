#pragma once

/*******************************************************************************
 * include/chaining.hpp
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
#include "utils/default_hash.hpp"
#include "utils/output.hpp"

#include "bucket.h"
#include "iterator_base.h"

namespace otm = utils_tm::out_tm;

namespace dysect
{
    struct triv_config
    { };

    template <class T>
    class iterator_incr;



    template<class K, class D, class HF = utils_tm::hash_tm::default_hash, class Conf = triv_config>
    class chaining
    {
    private:
        using this_type           = chaining<K,D,HF,Conf>;
        using hash_function_type  = HF;
        friend iterator_incr<this_type>;

    public:
        using size_type      = size_t;
        using key_type       = K;
        using mapped_type    = D;
        using iterator       = iterator_base<iterator_incr<this_type> >;
        using const_iterator = iterator_base<iterator_incr<this_type>, true>;

    private:
        using value_intern   = std::pair<key_type, mapped_type>;
        struct queue_item
        {
            queue_item*  next;
            value_intern element;
            queue_item(const std::pair<key_type, mapped_type>& element)
                : next(nullptr), element(element)
            { }
        };

    public:
        chaining(size_type cap, double alpha, size_type /*steps*/=0)
            : alpha(alpha), beta((alpha+1.)/2.), n(0),
              capacity((cap) ? cap*alpha : 2048*alpha),
              thresh  ((cap) ? cap*beta  : 2048*beta)
        {
            if (cap) table = std::make_unique<queue_item*[]>(capacity);
            std::fill(table.get(), table.get()+capacity, nullptr);
        }

        ~chaining()
        {
            if (!table) return;
            for (size_t i = 0; i<capacity; ++i)
            {
                auto curr = table[i];
                while (curr)
                {
                    auto temp = curr->next;
                    delete curr;
                    curr = temp;
                }
            }
        }

        chaining(const chaining&) = delete;
        chaining& operator=(const chaining&) = delete;

        chaining(chaining&& rhs)
            : alpha(rhs.alpha), beta(rhs.beta), n(rhs.n),
              capacity(rhs.capacity), thresh(rhs.thresh),
              hasher(rhs.hasher), table(std::move(rhs.table))
        {
            rhs.table = nullptr;
        }
        chaining& operator=(chaining&& rhs)
        {
            if (this == &rhs) return *this;
            this->~this_type();
            new (this) this_type(std::move(rhs));
            return *this;
        }

    private:
        /*** members that should become private at some point *********************/
        double     alpha;
        double     beta;
        size_type  n;
        size_type  capacity;
        size_type  thresh;
        hash_function_type  hasher;

        std::unique_ptr<queue_item*[]> table;

    public:
        // Basic Hash Table Functionality ******************************************
        iterator                  find  (const key_type& k);
        const_iterator            find  (const key_type& k) const;
        std::pair<iterator, bool> insert(const key_type& k, const mapped_type& d);
        std::pair<iterator, bool> insert(const value_intern& t);
        size_type                 erase (const key_type& k);
        int                 displacement(const key_type& k) const;
        void                      grow  ();

        // Easy use Accessors for std compliance ***********************************
        inline iterator           begin ();
        inline const_iterator     begin () const { return cbegin(); }
        inline const_iterator     cbegin() const;
        inline iterator           end   ()
        { return iterator      (nullptr, *this, nullptr, capacity); }
        inline const_iterator     end   () const { return cend(); }
        inline const_iterator     cend  () const
        { return const_iterator(nullptr, *this, nullptr, capacity); }

        mapped_type&              at    (const key_type& k);
        const mapped_type&        at    (const key_type& k) const;
        mapped_type&              operator[](const key_type& k);
        size_type                 count (const key_type& k) const;

        size_type                 get_capacity() const { return capacity; }

    private:
        void insert_intern(queue_item* item);

        // Easy iterators **********************************************************
        inline iterator       make_iterator (queue_item* item, size_type idx)
        { return iterator      (&(item->element), *this, item, idx); }
        inline const_iterator make_citerator(queue_item* item, size_type idx) const
        { return const_iterator(&(item->element), *this, item, idx); }

        // implementation specific functions (static polymorph) ********************
        inline size_type h(key_type k) const
        { return utils_tm::fastrange64(capacity, hasher(k)); }

        inline void inc_n()
        { if (++n > thresh) grow(); }
        inline void dec_n() { --n; }

        // Private helper function *************************************************

    public:
        inline static void print_init_header(otm::output_type& out)
        { out << otm::width(10) << "f_cap" << otm::width(16) << "overhead_bytes"; }
        inline void print_init_data  (otm::output_type& out)
        {
            long long overhead = (capacity + n)*sizeof(queue_item*)  // pointer
                + n*sizeof(std::pair<key_type, mapped_type>)    // elements
                - capacity*sizeof(std::pair<key_type, mapped_type>);
            out << otm::width(10) << capacity
                << otm::width(16) << overhead;
        }
    };



// Implementation of main functionality ****************************************

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::iterator
    chaining<K,D,H,C>::find(const key_type& k)
    {
        auto ind = h(k);

        auto curr = table[ind];

        while (curr)
        {
            if (curr->element.first == k) return make_iterator(curr, ind);

            curr = curr->next;
        }
        return end();
    }

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::const_iterator
    chaining<K,D,H,C>::find(const key_type& k) const
    {
        auto ind = h(k);

        auto curr = table[ind];

        while (curr)
        {
            if (curr->element.first == k) return make_citerator(curr, ind);

            curr = curr->next;
        }
        return cend();
    }

    template<class K, class D, class H, class C>
    inline std::pair<typename chaining<K,D,H,C>::iterator, bool>
    chaining<K,D,H,C>::insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    template<class K, class D, class H, class C>
    inline std::pair<typename chaining<K,D,H,C>::iterator, bool>
    chaining<K,D,H,C>::insert(const value_intern& t)
    {
        auto ind = h(t.first);

        auto curr = table[ind];
        if (!curr)
        {
            table[ind] = new queue_item(t);
            return std::make_pair(make_iterator(curr, ind), true);
        }
        while (true)
        {
            if (curr->element.first == t.first)
                return std::make_pair(make_iterator(curr, ind), false);
            if (!curr->next)
            {
                curr->next = new queue_item(t);
                return std::make_pair(make_iterator(curr, ind), true);
            }
            curr = curr->next;
        }
    }

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::size_type
    chaining<K,D,H,C>::erase(const key_type& k)
    {
        auto ind = h(k);

        auto curr = table[ind];
        if (!curr) return 0;
        if (curr->element.first == k)
        {
            table[ind] = curr->next;
            return 1;
        }

        auto prev = curr;
        curr = curr->next;

        while (curr)
        {
            if (curr->element.first == k)
            {
                prev->next = curr->next;
                return 1;
            }
            prev = curr;
            curr = curr->next;
        }
        return 0;
    }

    template<class K, class D, class H, class C>
    inline int
    chaining<K,D,H,C>::displacement(const key_type& k) const
    {
        auto ind = h(k);

        auto      curr  = table[ind];
        size_type steps = 0;

        while (curr)
        {
            if (curr->element.first == k) return steps;
            curr = curr->next;
            ++steps;
        }
        return -1;
    }


    template<class K, class D, class H, class C>
    inline void
    chaining<K,D,H,C>::grow()
    {
        this_type ntable(capacity, alpha);
        for (size_t i = 0; i < capacity; ++i)
        {
            while (table[i])
            {
                auto item = table[i];
                table[i] = item->next;
                ntable.insert_intern(item);
            }
        }
        (*this) = std::move(ntable);

    }

    template<class K, class D, class H, class C>
    inline void
    chaining<K,D,H,C>::insert_intern(queue_item* item)
    {
        auto ind   = h(item->element.first);
        item->next = table[ind];
        table[ind] = item;
    }

// Accessor Implementations ****************************************************

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::mapped_type&
    chaining<K,D,H,C>::at(const key_type& k)
    {
        auto a = find(k);
        if (a == end()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class K, class D, class H, class C>
    inline const typename chaining<K,D,H,C>::mapped_type&
    chaining<K,D,H,C>::at(const key_type& k) const
    {
        auto a = find(k);
        if (a == cend()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::mapped_type&
    chaining<K,D,H,C>::operator[](const key_type& k)
    {
        auto t = insert(k, mapped_type());
        return (*t.first).second;
    }

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::size_type
    chaining<K,D,H,C>::count(const key_type& k) const
    {
        return (find(k) != cend()) ? 1 : 0;
    }


    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::iterator
    chaining<K,D,H,C>::begin ()
    {
        for (size_t i = 0; i < capacity; ++i)
        {
            if (table[i])
                return make_iterator(table[i], i);
        }
        return end();
    }

    template<class K, class D, class H, class C>
    inline typename chaining<K,D,H,C>::const_iterator
    chaining<K,D,H,C>::cbegin() const
    {
        for (size_t i = 0; i < capacity; ++i)
        {
            if (table[i])
                return make_citerator(table[i], i);
        }
        return cend();
    }


// Iterator increment **********************************************************

    template<class Chaining>
    class iterator_incr
    {
    public:
        using table_type  = Chaining;

    private:
        using key_type    = typename table_type::key_type;
        using mapped_type = typename table_type::mapped_type;
        using ipointer    = std::pair<const key_type, mapped_type>*;
        using qitem       = typename table_type::queue_item;

    public:
        iterator_incr(const table_type& table, qitem* item, size_t index)
            : curr_item_ptr(item),
              table_ptr(&table.table[index]),
              end_table_ptr(&table.table[table.capacity-1])
        { }
        iterator_incr(const iterator_incr&) = default;
        iterator_incr& operator=(const iterator_incr&) = default;

        ipointer next(ipointer cur)
        {
            if (curr_item_ptr->next)
            {
                curr_item_ptr = curr_item_ptr->next;
                return &(curr_item_ptr->element);
            }
            while (table_ptr < end_table_ptr)
            {
                ++table_ptr;
                if (*table_ptr)
                {
                    curr_item_ptr = *table_ptr;
                    return &(curr_item_ptr->element);
                }
            }
            return nullptr;
        }

    private:
        qitem*  curr_item_ptr;
        qitem** table_ptr;
        qitem** end_table_ptr;
    };

}
