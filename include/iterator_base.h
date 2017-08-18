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


// template <class Table, class Ptr>
// class iterator_incr
// {
// private:
//     using Pointer_t = Ptr;
// public:
//     iterator_incr(const Table&) { }

//     iterator_incr(const iterator_incr&) = default;
//     iterator_incr& operator=(const iterator_incr&) = default;

//     Pointer_t next(Pointer_t) { return nullptr; }
// };

template <class Increment, bool is_const = false>
class IteratorBase
{
private:
    using Table_t      = typename Increment::Table_t;

    using key_type     = typename Table_t::key_type;
    using mapped_type  = typename Table_t::mapped_type;
    using value_intern = std::pair<      key_type, mapped_type>;
    using value_table  = std::pair<const key_type, mapped_type>;
    using cval_intern  = typename std::conditional<is_const, const value_intern , value_intern >::type;

public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::conditional<is_const, const value_table, value_table>::type;
    using reference  = value_type&;
    using pointer    = value_type*;
    using iterator_category = std::forward_iterator_tag;
    using Incr_t   = Increment;

    template<class T, bool b>
    friend void swap(IteratorBase<T,b>& l, IteratorBase<T,b>& r);
    template<class T, bool b>
    friend bool operator==(const IteratorBase<T,b>& l, const IteratorBase<T,b>& r);
    template<class T, bool b>
    friend bool operator!=(const IteratorBase<T,b>& l, const IteratorBase<T,b>& r);


    // Constructors ************************************************************

    IteratorBase(cval_intern* pair_, const Table_t& table)
        : ptr(reinterpret_cast<pointer>(pair_)), incr(table) { }

    IteratorBase(const IteratorBase& rhs) : ptr(rhs.ptr), incr(rhs.incr) { }
    IteratorBase& operator=(const IteratorBase& r)
    { ptr = r.ptr; incr = r.incr; return *this; }

    ~IteratorBase() = default;


    // Basic Iterator Functionality

    IteratorBase& operator++(int = 0) { ptr = incr.next(ptr); return *this; }
    reference operator* () const { return *ptr; }
    pointer   operator->() const { return  ptr; }

    bool operator==(const IteratorBase& rhs) const { return ptr == rhs.ptr; }
    bool operator!=(const IteratorBase& rhs) const { return ptr != rhs.ptr; }

private:
    pointer ptr;
    Incr_t  incr;
};
