#pragma once

#include "cuckoo_base.h"

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooSimple : public CuckooTraits<CuckooSimple<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooSimple<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using size_type      = typename CuckooTraits<This_t>::size_type;
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

private:
    using value_intern = std::pair<key_type, mapped_type>;

public:
    CuckooSimple(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          beta((size_constraint + 1.)/2.)
    {
        n_buckets   = size_type(double(cap)*size_constraint)/bs;
        n_buckets   = std::max<size_type>(n_buckets, 256);

        capacity    = n_buckets*bs;
        grow_thresh = beta*std::max<size_type>(256ull, cap);

        factor      = double(n_buckets)/double(1ull<<32);

        table       = std::make_unique<Bucket_t[]>(n_buckets);
    }

    CuckooSimple(const CuckooSimple&) = delete;
    CuckooSimple& operator=(const CuckooSimple&) = delete;

    CuckooSimple(CuckooSimple&&) = default;
    CuckooSimple& operator=(CuckooSimple&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    static constexpr size_type bs = CuckooTraits<This_t>::bs;
    static constexpr size_type nh = CuckooTraits<This_t>::nh;
    static constexpr size_type tl = 1;

    size_type n_buckets;
    double beta;
    double factor;

    std::unique_ptr<Bucket_t[]> table;
    std::vector<value_intern> grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
                     : std::make_pair(0,nullptr);
    }

    inline iterator begin()
    {
        auto temp = make_iterator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    inline const_iterator cbegin() const
    {
        auto temp = make_citerator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    // Functions for finding buckets *******************************************

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_type i = 0; i < nh; ++i)
        {
            mem[i] = getBucket(h,i);
        }
    }

    inline Bucket_t* getBucket(Hashed_t h, size_type i) const
    {
        size_type l = Ext::loc(h,i) * factor;
        return &(table[l]);
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        if (grow_buffer.size()) return;
        size_type nsize = size_type(double(n)*alpha) / bs;
        nsize        = std::max(nsize, n_buckets+1);
        capacity     = nsize*bs;
        double nfactor = double(nsize)/double(1ull << 32);
        size_type nthresh = n * beta;

        //std::cout << n << " " << n_buckets << " -> " << nsize << std::endl;

        auto ntable = std::make_unique<Bucket_t[]>(nsize);
        migrate(ntable, nfactor);

        n_buckets   = nsize;
        table       = std::move(ntable);
        grow_thresh = nthresh;
        factor      = nfactor;
        if (grow_buffer.size()) finalize_grow();
    }

    inline void migrate(std::unique_ptr<Bucket_t[]>& target, double nfactor)
    {
        for (size_type i = 0; i < n_buckets; ++i)
        {
            Bucket_t& curr = table[i];

            for (size_type j = 0; j < bs; ++j)
            {
                auto e = curr.elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);
                for (size_type ti = 0; ti < nh; ++ti)
                {
                    if (i == size_type(Ext::loc(hash, ti)*factor))
                    {
                        if (! target[Ext::loc(hash, ti) * nfactor].insert(e.first, e.second))
                        {
                            grow_buffer.push_back(e);
                        }
                        break;
                    }
                }
            }
        }
    }

    inline void finalize_grow()
    {
        size_type temp = n;
        for (auto& e : grow_buffer)
        {
            insert(e);
        }
        n = temp;
        std::vector<value_intern> tbuf;
        //grow_buffer.clear();
        std::swap(tbuf, grow_buffer);
    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooSimple<K,D,HF,Conf> >
{
public:
    using Specialized_t = CuckooSimple<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Config_t      = Conf;

    using size_type     = size_t;
    using key_type      = K;
    using mapped_type   = D;

    static constexpr size_type tl = 1;
    static constexpr size_type bs = Conf::bs;
    static constexpr size_type nh = Conf::nh;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooSimple<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooSimple<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using pointer   = std::pair<const K,D>*;
    static constexpr size_type bs = Conf::bs;

public:
    iterator_incr(const Table_t& table_)
        : end_ptr(reinterpret_cast<pointer>
                  (&table_.table[table_.n_buckets-1].elements[bs-1]))
    { }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    pointer next(pointer cur)
    {
        while (cur < end_ptr)
        {
            if ((++cur)->first) return cur;
        }
        return nullptr;
    }

private:
    pointer end_ptr;
};












// *****************************************************************************
// IN PLACE GROWING ************************************************************
// *****************************************************************************

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooInPlace : public CuckooTraits<CuckooInPlace<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooInPlace<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using size_type      = typename CuckooTraits<This_t>::size_type;
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

private:
    using value_intern = std::pair<key_type, mapped_type>;


    static constexpr size_type bs = CuckooTraits<This_t>::bs;
    static constexpr size_type nh = CuckooTraits<This_t>::nh;
    static constexpr size_type tl = 1;

    static constexpr size_type max_size     = 10ull << 30;

public:
    CuckooInPlace(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          beta((size_constraint + 1.)/2.)
    {
        auto temp = reinterpret_cast<Bucket_t*>(operator new (max_size));
        table = std::unique_ptr<Bucket_t[]>(temp);

        n_buckets   = size_type(double(cap)*size_constraint)/bs;
        n_buckets   = std::max<size_type>(n_buckets, 256);

        capacity    = n_buckets*bs;
        grow_thresh = beta*std::max<size_type>(256ull, cap);

        factor      = double(n_buckets)/double(1ull<<32);

        //table       = std::make_unique<Bucket_t[]>(n_buckets);
        std::fill(table.get(), table.get()+n_buckets, Bucket_t());
    }

    CuckooInPlace(const CuckooInPlace&) = delete;
    CuckooInPlace& operator=(const CuckooInPlace&) = delete;

    CuckooInPlace(CuckooInPlace&&) = default;
    CuckooInPlace& operator=(CuckooInPlace&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    size_type n_buckets;
    double    beta;
    double    factor;

    std::unique_ptr<Bucket_t[]> table;
    std::vector<value_intern>   grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (! i) ? std::make_pair(n_buckets, table.get())
                     : std::make_pair(0,nullptr);
    }

    inline iterator begin()
    {
        auto temp = make_iterator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    inline const_iterator cbegin() const
    {
        auto temp = make_citerator(&table[0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    // Functions for finding buckets *******************************************

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_type i = 0; i < nh; ++i)
        {
            mem[i] = getBucket(h,i);
        }
    }

    inline Bucket_t* getBucket(Hashed_t h, size_type i) const
    {
        size_type l = Ext::loc(h,i) * factor;
        return &(table[l]);
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        // std::cout << "!" << std::flush;

        if (grow_buffer.size()) return;
        size_type nsize   = size_type(double(n)*alpha) / bs;
        nsize             = std::max(nsize, n_buckets+1);
        capacity          = nsize*bs;
        double    nfactor = double(nsize)/double(1ull << 32);
        size_type nthresh = n * beta;

        //std::cout << n << " " << n_buckets << " -> " << nsize << std::endl;

        //auto ntable = std::make_unique<Bucket_t[]>(nsize);
        std::fill(table.get()+n_buckets, table.get()+nsize, Bucket_t());

        migrate(nfactor);

        n_buckets   = nsize;
        //table       = std::move(ntable);
        grow_thresh = nthresh;
        factor      = nfactor;
        if (grow_buffer.size()) finalize_grow();

        // std::cout << n_buckets << " " << std::flush;
        // std::cout << grow_thresh << "!" << std::endl;
    }

    inline void migrate(double nfactor)
    {
        for (int i = n_buckets; i >= 0; --i)
        {
            Bucket_t& curr = table[i];

            for (size_type j = 0; j < bs; ++j)
            {
                auto e = curr.elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_type ti = 0; ti < nh; ++ti)
                {
                    if (i == int (Ext::loc(hash, ti)*factor))
                    {
                        auto nbucket = Ext::loc(hash,ti) * nfactor;
                        if (    (i == nbucket)
                             || (! table[nbucket].insert(e.first, e.second)) )
                        {
                            grow_buffer.push_back(e);
                        }
                        break;
                    }
                }
            }
            curr = Bucket_t();
        }
    }

    inline void finalize_grow()
    {
        size_type temp = n;
        for (auto& e : grow_buffer)
        {
            insert(e);
        }
        n = temp;
        std::vector<value_intern> tbuf;
        //grow_buffer.clear();
        std::swap(tbuf, grow_buffer);
    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooInPlace<K,D,HF,Conf> >
{
public:
    using Specialized_t = CuckooInPlace<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Config_t      = Conf;

    using size_type     = size_t;
    using key_type      = K;
    using mapped_type   = D;

    static constexpr size_type tl = 1;
    static constexpr size_type bs = Conf::bs;
    static constexpr size_type nh = Conf::nh;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooInPlace<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooInPlace<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using pointer   = std::pair<const K,D>*;
    static constexpr size_type bs = Conf::bs;

public:
    iterator_incr(const Table_t& table_)
        : end_ptr(reinterpret_cast<pointer>
                  (&table_.table[table_.n_buckets-1].elements[bs-1]))
    { }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    pointer next(pointer cur)
    {
        while (cur < end_ptr)
        {
            if ((++cur)->first) return cur;
        }
        return nullptr;
    }

private:
    pointer end_ptr;
};
