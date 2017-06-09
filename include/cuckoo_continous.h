#pragma once

#include "cuckoo_base.h"
#include "cbucket.h"

template<class K0, class D0, class HF0, class Conf0>
class CuckooIndependentBase;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class cuckoo_overlap : public CuckooTraits<cuckoo_overlap<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = cuckoo_overlap<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

    using indep_var = CuckooIndependentBase<K,D,HF,Conf>;
    friend indep_var;
public:
    using size_type      = typename CuckooTraits<This_t>::size_type;
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;

private:
    using value_intern = std::pair<key_type, mapped_type>;

public:
    cuckoo_overlap(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          beta((size_constraint + 1.)/2.)
    {
        n_subbuckets   = size_type(double(cap)*size_constraint)/sbs;
        n_subbuckets   = std::max<size_type>(n_subbuckets, 256);

        capacity    = n_subbuckets*sbs;
        grow_thresh = beta*std::max<size_type>(256ull, cap);

        factor      = double(n_subbuckets+1-(bs/sbs))/double(1ull<<32);

        table       = std::make_unique<value_type[]>(n_subbuckets*sbs);
    }

    cuckoo_overlap(const cuckoo_overlap&) = delete;
    cuckoo_overlap& operator=(const cuckoo_overlap&) = delete;

    cuckoo_overlap(cuckoo_overlap&&) = default;
    cuckoo_overlap& operator=(cuckoo_overlap&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    static constexpr size_type sbs = CuckooTraits<This_t>::sbs;
    static constexpr size_type bs  = CuckooTraits<This_t>::bs;
    static constexpr size_type nh  = CuckooTraits<This_t>::nh;
    static constexpr size_type tl  = 1;

    size_type n_subbuckets;
    double beta;
    double factor;

    std::unique_ptr<value_intern[]> table;
    std::vector<value_intern> grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

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
        return reinterpret_cast<Bucket_t*>(&(table[l*sbs]));
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        if (grow_buffer.size()) return;
        size_type nsize = size_type(double(n)*alpha) / sbs;
        nsize           = std::max(nsize, n_subbuckets+1);
        capacity        = nsize*bs;
        double nfactor  = double(nsize+1-(bs/sbs))/double(1ull << 32);
        size_type nthresh = n * beta;

        //std::cout << n << " " << n_buckets << " -> " << nsize << std::endl;

        auto ntable     = std::make_unique<value_intern[]>(nsize*sbs);
        migrate(ntable, nfactor);

        n_subbuckets    = nsize;
        table           = std::move(ntable);
        grow_thresh     = nthresh;
        factor          = nfactor;
        if (grow_buffer.size()) finalize_grow();
    }

    inline void migrate(std::unique_ptr<Bucket_t[]>& target, double nfactor)
    {
        for (size_type i = 0; i < capacity; ++i)
        {
            auto e = table[i];
            if (! e.first) continue;
            auto hash = hasher(e.first);

            for (size_type ti = 0; ti < nh; ++ti)
            {
                size_type off_ti = size_type(Ext::loc(hash,ti) * factor)*sbs;
                if (i >= off_ti && i < off_ti+bs)
                {
                    if (! reinterpret_cast<Bucket_t*>
                        (&target[size_type(Ext::loc(hash, ti)*nfactor)*sbs])
                        .insert(e))
                    {
                        grow_buffer.push_back(e);
                    }
                }
            }

            // Bucket_t& curr = table[i];

            // for (size_type j = 0; j < bs; ++j)
            // {
            //     auto e = curr.elements[j];
            //     if (! e.first) break;
            //     auto hash = hasher(e.first);
            //     for (size_type ti = 0; ti < nh; ++ti)
            //     {
            //         if (i == size_type(Ext::loc(hash, ti)*factor))
            //         {
            //             if (! target[Ext::loc(hash, ti) * nfactor].insert(e.first, e.second))
            //             {
            //                 grow_buffer.push_back(e);
            //             }
            //             break;
            //         }
            //     }
            // }
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
class CuckooTraits<cuckoo_overlap<K,D,HF,Conf> >
{
public:
    using Specialized_t = cuckoo_overlap<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Config_t      = Conf;

    using size_type     = size_t;
    using key_type      = K;
    using mapped_type   = D;

    static constexpr size_type tl  = 1;
    static constexpr size_type bs  = Conf::bs;
    static constexpr size_type sbs = Conf::sbs;
    static constexpr size_type nh  = Conf::nh;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = CBucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<cuckoo_overlap<K,D,HF,Conf> >
{
public:
    using Table_t   = cuckoo_overlap<K,D,HF,Conf>;
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
class cuckoo_inplace_overlap : public CuckooTraits<cuckoo_inplace_overlap<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = cuckoo_inplace_overlap<K,D,HF,Conf>;
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

    static constexpr size_type bs  = CuckooTraits<This_t>::bs;
    static constexpr size_type sbs = CuckooTraits<This_t>::sbs;
    static constexpr size_type nh  = CuckooTraits<This_t>::nh;
    static constexpr size_type tl  = 1;

    static constexpr size_type max_size     = 10ull << 30;

public:
    cuckoo_inplace_overlap(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          beta((size_constraint + 1.)/2.)
    {
        auto temp    = reinterpret_cast<value_intern*>(operator new (max_size));
        table        = std::unique_ptr<value_intern[]>(temp);

        n_subbuckets = size_type(double(cap)*size_constraint)/sbs;
        n_subbuckets = std::max<size_type>(n_subbuckets, 256);

        capacity     = n_subbuckets*sbs;
        grow_thresh  = beta*std::max<size_type>(256ull, cap);

        factor       = double(n_subuckets+1-(bs/sbs))/double(1ull<<32);

        std::fill(table.get(), table.get()+capacity, value_intern());
    }

    cuckoo_inplace_overlap(const cuckoo_inplace_overlap&) = delete;
    cuckoo_inplace_overlap& operator=(const cuckoo_inplace_overlap&) = delete;

    cuckoo_inplace_overlap(cuckoo_inplace_overlap&&) = default;
    cuckoo_inplace_overlap& operator=(cuckoo_inplace_overlap&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    size_type n_subbuckets;
    double    beta;
    double    factor;

    std::unique_ptr<value_intern[]> table;
    std::vector<value_intern>   grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

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
        return reinterpret_cast<Bucket_t*>(&(table[l*sbs]));
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        if (grow_buffer.size()) return;
        size_type nsize   = size_type(double(n)*alpha) / sbs;
        nsize             = std::max(nsize, n_subbuckets+1);
        capacity          = nsize*sbs;
        double    nfactor = double(nsize+1-(bs/sbs))/double(1ull << 32);
        size_type nthresh = n * beta;

        std::fill(table.get()+capacity, table.get()+(nsize*sbs), value_intern());

        migrate(nfactor);

        n_subbuckets      = nsize;
        grow_thresh       = nthresh;
        factor            = nfactor;
        if (grow_buffer.size()) finalize_grow();
    }

    inline void migrate(double nfactor)
    {
        for (int i = capacity-1; i >= 0; --i)
        {
            auto e = table[i];
            if  (!e.first) continue;
            table[i] = value_intern();
            auto hash = hasher(e.first);

            for (size_type ti = 0; ti < nh; ++ti)
            {
                size_type off_i = size_type(Ext::loc(hash, ti)*factor)*sbs;
                if (i >= int(off_i) && i < int(off_i + bs))
                {
                    size_type nbucket = size_type(Ext::loc(hash,ti) * nfactor)*sbs;
                    if (nbucket >= i ||
                        ! reinterpret_cast<Bucket_t*>(&table[nbucket])->insert(e));
                    {
                        grow_buffer.push_back(e);
                    }
                    break;
                }
            }

            // Bucket_t& curr = table[i];

            // for (size_type j = 0; j < bs; ++j)
            // {
            //     auto e = curr.elements[j];
            //     if (! e.first) break;
            //     auto hash = hasher(e.first);

            //     for (size_type ti = 0; ti < nh; ++ti)
            //     {
            //         if (i == int (Ext::loc(hash, ti)*factor))
            //         {
            //             auto nbucket = Ext::loc(hash,ti) * nfactor;
            //             if (    (i == nbucket)
            //                  || (! table[nbucket].insert(e.first, e.second)) )
            //             {
            //                 grow_buffer.push_back(e);
            //             }
            //             break;
            //         }
            //     }
            // }
            // curr = Bucket_t();
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
class CuckooTraits<cuckoo_inplace_overlap<K,D,HF,Conf> >
{
public:
    using Specialized_t = cuckoo_inplace_overlap<K,D,HF,Conf>;
    using Base_t        = CuckooMultiBase<Specialized_t>;
    using Config_t      = Conf;

    using size_type     = size_t;
    using key_type      = K;
    using mapped_type   = D;

    static constexpr size_type tl = 1;
    static constexpr size_type bs = Conf::bs;
    static constexpr size_type nh = Conf::nh;
    static constexpr size_type sbs= Conf::sbs;

    using Hasher_t      = Hasher<K, HF, 0, nh, true, true>;
    using Bucket_t      = Bucket<K,D,bs>;

};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<cuckoo_inplace_overlap<K,D,HF,Conf> >
{
public:
    using Table_t   = cuckoo_inplace_overlap<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using pointer   = std::pair<const K,D>*;

public:
    iterator_incr(const Table_t& table_)
        : end_ptr(reinterpret_cast<pointer>
                  (&table_.table[table_.capacity-1]))
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
