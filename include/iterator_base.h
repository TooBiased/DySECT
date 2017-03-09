
template <class Table>
class IteratorBase
    {
        using Table_t = Table;
        using Key     = typename Table_t::Key;
        using Data    = typename Table_t::Data;
        using Pair_t  = std::pair<Key, Data>;
    public:
        using value_type = Pair_t;
        using reference  = value_type&;
        using pointer    = value_type*;
        using iterator_category = std::forward_iterator_tag;

        template<class T>
        friend void swap(IteratorBase<T>& lhs, IteratorBase<T>& rhs);
        template<class T>
        friend bool operator==(const IteratorBase<T>& lhs, const IteratorBase<T>& rhs);
        template<class T>
        friend bool operator!=(const IteratorBase<T>& lhs, const IteratorBase<T>& rhs);

        IteratorBase(Pair_t* pair_) : pair(pair_) { }
        IteratorBase(const IteratorBase& rhs)   : pair(rhs.pair) { }
        IteratorBase& operator=(const IteratorBase& rhs)
        { pair = rhs.pair; return *this; }
        ~IteratorBase() { }

        IteratorBase& operator++(int = 0) { return end(); }
        reference operator* () const  { return *pair; }
        pointer   operator->() const  { return pair; }

        bool operator==(const IteratorBase& rhs) { return pair == rhs.pair; }
        bool operator!=(const IteratorBase& rhs) { return pair != rhs.pair; }

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
        Pair_t* pair;
        static constexpr IteratorBase end() { return IteratorBase(nullptr); }
    };
