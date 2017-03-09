#pragma once

template <class Table, bool is_const = false>
class IteratorBase
{
private:
    using Table_t = Table;
    using Key     = typename Table_t::Key;
    using Data    = typename Table_t::Data;
    using stored  = std::pair<      Key, Data>;
    using pair    = std::pair<const Key, Data>;
    using c_stored = typename std::conditional<is_const, const stored, stored>::type;

public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::conditional<is_const, const pair, pair>::type;
    using reference  = value_type&;
    using pointer    = value_type*;
    using iterator_category = std::forward_iterator_tag;

    template<class T, bool b>
    friend void swap(IteratorBase<T,b>& l, IteratorBase<T,b>& r);
    template<class T, bool b>
    friend bool operator==(const IteratorBase<T,b>& l, const IteratorBase<T,b>& r);
    template<class T, bool b>
    friend bool operator!=(const IteratorBase<T,b>& l, const IteratorBase<T,b>& r);

    //IteratorBase(pointer pair_ = nullptr)
    //    : ptr(pair_) { }
    IteratorBase(c_stored* pair_)
        : ptr(reinterpret_cast<pointer>(pair_)) { }

    IteratorBase(const IteratorBase& rhs)   : ptr(rhs.ptr) { }
    IteratorBase& operator=(const IteratorBase& r) { ptr = r.ptr; return *this; }

    ~IteratorBase() = default;

    IteratorBase& operator++(int = 0) { return end(); }
    reference operator* () const { return *ptr; }
    pointer   operator->() const { return  ptr; }

    bool operator==(const IteratorBase& rhs) const { return ptr == rhs.ptr; }
    bool operator!=(const IteratorBase& rhs) const { return ptr != rhs.ptr; }

    // bool valid()
    // {
    //     if (!pair)                    return false;
    //     if (table.check_verison(ver)) return false;
    //     if (pair->key != key)         return false;
    //     return true;
    // }
    // bool restore()
    // {
    //     *this = table.findIt(key);
    // }
    // table_t table&;
    // size_t  ver;
    // key_t   key;
    pointer ptr;

    static constexpr IteratorBase end() { return IteratorBase(nullptr); }
};
