#pragma once

/*******************************************************************************
 * include/iterator_base.h
 *
 * iterator_base is a basic iterator implementation, that is
 * independent from the used table. Appart from the increment
 * operator, which uses a table specific increment function.
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

    template <class Increment, bool is_const = false>
    class iterator_base
    {
    private:
        using table_type      = typename Increment::table_type;

        using key_type     = typename table_type::key_type;
        using mapped_type  = typename table_type::mapped_type;
        using value_intern = std::pair<      key_type, mapped_type>;
        using value_table  = std::pair<const key_type, mapped_type>;
        using cval_intern  = typename std::conditional<is_const, const value_intern , value_intern >::type;

    public:
        using difference_type = std::ptrdiff_t;
        using value_type = typename std::conditional<is_const, const value_table, value_table>::type;
        using reference  = value_type&;
        using pointer    = value_type*;
        using iterator_category = std::forward_iterator_tag;
        using increment_type   = Increment;

        template<class T, bool b>
        friend void swap(iterator_base<T,b>& l, iterator_base<T,b>& r);
        template<class T, bool b>
        friend bool operator==(const iterator_base<T,b>& l, const iterator_base<T,b>& r);
        template<class T, bool b>
        friend bool operator!=(const iterator_base<T,b>& l, const iterator_base<T,b>& r);


        // Constructors ************************************************************
        template<class ... Args>
        iterator_base(cval_intern* pair_, Args&& ... args)
            : ptr(reinterpret_cast<pointer>(pair_)), incr(args...) { }

        iterator_base(const iterator_base& rhs) : ptr(rhs.ptr), incr(rhs.incr) { }
        iterator_base& operator=(const iterator_base& r)
        { ptr = r.ptr; incr = r.incr; return *this; }

        ~iterator_base() = default;


        // Basic Iterator Functionality

        iterator_base& operator++(int) { ptr = incr.next(ptr); return *this; }
        reference operator* () const { return *ptr; }
        pointer   operator->() const { return  ptr; }

        bool operator==(const iterator_base& rhs) const { return ptr == rhs.ptr; }
        bool operator!=(const iterator_base& rhs) const { return ptr != rhs.ptr; }

    private:
        pointer        ptr;
        increment_type incr;
    };

}
