#pragma once

#include <cmath>
#include "cuckoo_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooHomogeneous2L : public CuckooTraits<CuckooHomogeneous2L<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooHomogeneous2L<K,D,HF,Conf>;
    using Base_t         = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t       = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t       = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t       = typename Hasher_t::Hashed_t;
    using Ext            = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using key_type       = typename CuckooTraits<This_t>::key_type;
    using mapped_type    = typename CuckooTraits<This_t>::mapped_type;
    using iterator       = typename Base_t::iterator;
    using const_iterator = typename Base_t::const_iterator;


    CuckooHomogeneous2L(size_t cap = 0      , double size_constraint = 1.1,
                       size_t dis_steps = 0, size_t seed = 0)
        : Base_t(size_constraint, dis_steps, seed),
          ll_size(std::floor(double(cap) * size_constraint / double(tl*bs))),
          beta((size_constraint+1.)/2.),
          factor(double(ll_size)/double(1ull << (32-ct_log(tl))))
    {
        for (size_t i = 0; i < tl; ++i)
        {
            ll_table[i] = std::make_unique<Bucket_t[]>(ll_size);
        }
        capacity    = tl * ll_size * bs;
        grow_thresh = double(capacity) / beta;
    }

    CuckooHomogeneous2L(const CuckooHomogeneous2L&) = delete;
    CuckooHomogeneous2L& operator=(const CuckooHomogeneous2L&) = delete;

    CuckooHomogeneous2L(CuckooHomogeneous2L&& rhs)
        : Base_t(std::move(rhs)), ll_size(rhs.ll_size), beta(rhs.beta),
          factor(rhs.factor)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            ll_table[i] = std::move(rhs.ll_table[i]);
        }
    }

    CuckooHomogeneous2L& operator=(CuckooHomogeneous2L&& rhs)
    {
        Base_t::operator=(std::move(rhs));

        ll_size= rhs.ll_size;
        beta   = rhs.beta;
        factor = rhs.factor;

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(ll_table[i], rhs.ll_table[i]);
        }
        return *this;
    }

    ~CuckooHomogeneous2L() = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;

    size_t                      ll_size;
    double                      beta;
    double                      factor;

    std::unique_ptr<Bucket_t[]> ll_table[tl];
    std::vector<std::pair<key_type, mapped_type> > grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(ll_size, ll_table[i].get())
                        : std::make_pair(0,nullptr);
    }

    inline iterator begin()
    {
        auto temp = make_iterator(&ll_table[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    inline const_iterator cbegin() const
    {
        auto temp = make_citerator(&ll_table[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }


private:
    // Functions for finding buckets *******************************************

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h,i);
    }

    inline Bucket_t* getBucket(Hashed_t h, size_t i) const
    {
        size_t tab = Ext::tab(h,i);
        return &(ll_table[tab][Ext::loc(h,i) * factor]);
    }



    // Size changes (GROWING) **************************************************

    void grow()
    {
        size_t nll_size = std::floor(double(n)*alpha / double(tl*bs));
        nll_size = std::max(nll_size, ll_size+1);
        double nfactor  = double(nll_size)/double(1ull << (32-ct_log(tl)));

        for (size_t i = 0; i < tl; ++i)
        {
            auto ntable = std::make_unique<Bucket_t[]>(nll_size);
            migrate(i, ntable, nfactor);
            ll_table[i] = std::move(ntable);
        }

        ll_size  = nll_size;
        factor   = nfactor;
        capacity = ll_size*tl*bs;
        grow_thresh   = double(n)*beta;
        if (grow_buffer.size()) finalize_grow();
    }

    inline void migrate(size_t ind, std::unique_ptr<Bucket_t[]>& target, double tfactor)
    {
        //auto& source = ll_table[ind];

        for(size_t i = 0; i < ll_size; ++i)
        {
            Bucket_t& bucket = ll_table[ind][i];

            for (size_t j = 0; j < bs; ++j)
            {
                auto element = bucket.elements[j];
                if (! element.first) break;
                auto hash = hasher(element.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if ((ind == Ext::tab(hash,ti)) &&
                        (i   == size_t(Ext::loc(hash, ti)*factor)))
                    {
                        if (!target[Ext::loc(hash,ti) * tfactor].insert(element))
                            grow_buffer.push_back(element);
                    }
                }
            }
        }
    }

    inline void finalize_grow()
    {
        size_t temp = n;
        for (auto& e : grow_buffer)
        {
            insert(e);
        }
        n = temp;
        grow_buffer.clear();
    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooHomogeneous2L<K,D,HF,Conf> >
{
public:
    using Specialized_t  = CuckooHomogeneous2L<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Config_t       = Conf;

    using key_type       = K;
    using mapped_type    = D;

    static constexpr size_t tl = Config_t::tl;
    static constexpr size_t bs = Config_t::bs;
    static constexpr size_t nh = Config_t::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;
};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooHomogeneous2L<K,D,HF,Conf> >
{
private:
    using ipointer = std::pair<const K,D>*;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t bs = Conf::bs;

public:
    using Table_t   = CuckooHomogeneous2L<K,D,HF,Conf>;

    iterator_incr(const Table_t& table_)
        : table(table_), end_tab(nullptr), tab(tl + 1)
    { }
    iterator_incr(const iterator_incr&) = default;
    iterator_incr& operator=(const iterator_incr&) = default;

    ipointer next(ipointer cur)
    {
        if (tab > tl) initialize_tab(cur);

        auto temp = cur+1;
        if (temp > end_tab)
        { temp = overflow_tab(); if (!temp) return nullptr; }

        while (!temp->first)
        {
            if (++temp > end_tab) return overflow_tab();
        }
        return temp;
    }

private:
    const Table_t& table;
    ipointer       end_tab;
    size_t         tab;

    ipointer overflow_tab()
    {
        if (++tab >= tl) return nullptr;
        size_t size = table.ll_size - 1;
        end_tab = reinterpret_cast<ipointer>(&table.ll_table[tab][size].elements[bs-1]);
        return reinterpret_cast<ipointer>(&table.ll_table[tab][0].elements[0]);
    }

    void initialize_tab(ipointer ptr)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            ipointer tab_b_ptr = reinterpret_cast<ipointer>
                (&table.ll_table[i][0].elements[0]);
            ipointer tab_e_ptr = tab_b_ptr + table.ll_size*bs;

            if (tab_b_ptr <= ptr && ptr < tab_e_ptr)
            {
                tab = i;
                end_tab = tab_e_ptr;
                return;
            }
        }
    }

};
