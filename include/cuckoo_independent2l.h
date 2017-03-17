#pragma once

#include <cmath>
#include "cuckoo_base.h"

template<class T>
class CuckooTraits;

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooIndependent2L
    : public CuckooTraits<CuckooIndependent2L<K,D,HF,Conf> >::Base_t
{
private:
    using This_t             = CuckooIndependent2L<K,D,HF,Conf>;
    using Base_t             = typename CuckooTraits<This_t>::Base_t;
    using Bucket_t           = typename CuckooTraits<This_t>::Bucket_t;
    using Hasher_t           = typename CuckooTraits<This_t>::Hasher_t;
    using Hashed_t           = typename Hasher_t::Hashed_t;
    using Ext                = typename Hasher_t::Extractor_t;

    friend Base_t;
    friend iterator_incr<This_t>;

public:
    using key_type           = typename CuckooTraits<This_t>::key_type;
    using mapped_type        = typename CuckooTraits<This_t>::mapped_type;
    using iterator           = typename Base_t::iterator;
    using const_iterator     = typename Base_t::const_iterator;
    using insert_return_type = typename Base_t::insert_return_type;

private:
    using value_intern   = std::pair<key_type, mapped_type>;

public:
    CuckooIndependent2L(size_t cap = 0      , double size_constraint = 1.1,
                        size_t dis_steps = 0, size_t seed = 0)
        : Base_t(size_constraint, dis_steps, seed), beta((1.0+size_constraint)/2.0)
    {
        size_t lsize  = std::floor(cap * size_constraint / double(tl*bs));
        lsize  = std::max(lsize, 256ul);
        double factor = double(lsize)/fac_div;
        size_t grow_thresh = double(lsize) / beta;

        for (size_t i = 0; i < tl; ++i)
        {
            ll_tab[i]    = std::make_unique<Bucket_t[]>(lsize);
            ll_size[i]   = lsize;
            ll_elem[i]   = 0;
            ll_thresh[i] = grow_thresh;
            ll_factor[i] = factor;
        }

        capacity    = bs * tl * lsize;
    }

    CuckooIndependent2L(const CuckooIndependent2L&) = delete;
    CuckooIndependent2L& operator=(const CuckooIndependent2L&) = delete;

    CuckooIndependent2L(CuckooIndependent2L&& rhs)
        : Base_t(std::move(rhs)), beta(rhs.beta)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            ll_tab[i]    = std::move(rhs.ll_tab[i]);
            ll_size[i]   = rhs.ll_size[i];
            ll_elem[i]   = rhs.ll_elem[i];
            ll_thresh[i] = rhs.ll_thresh[i];
            ll_factor[i] = rhs.ll_factor[i];
        }
    }

    CuckooIndependent2L& operator=(CuckooIndependent2L&& rhs)
    {
        Base_t::operator=(std::move(rhs));
        beta            = rhs.beta;

        for (size_t i = 0; i < tl; ++i)
        {
            std::swap(ll_tab[i]   , rhs.ll_tab[i]);
            std::swap(ll_size[i]  , rhs.ll_size[i]);
            std::swap(ll_elem[i]  , rhs.ll_elem[i]);
            std::swap(ll_thresh[i], rhs.ll_thresh[i]);
            std::swap(ll_factor[i], rhs.ll_factor[i]);
        }
        return *this;
    }

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::alpha;
    using Base_t::hasher;

    static constexpr size_t bs = CuckooTraits<This_t>::bs;
    static constexpr size_t tl = CuckooTraits<This_t>::tl;
    static constexpr size_t nh = CuckooTraits<This_t>::nh;
    static constexpr double fac_div = double (1ull << (32 - ct_log(tl)));

    double beta;
    size_t     ll_size  [tl];
    size_t     ll_elem  [tl];
    size_t     ll_thresh[tl];
    double     ll_factor[tl];
    std::unique_ptr<Bucket_t[]> ll_tab[tl];
    std::vector<value_intern> grow_buffer;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    // Specialized Funcitions (to keep per table counts) ***********************

    inline insert_return_type insert(const key_type k, const mapped_type d)
    {
        return insert(std::make_pair(k,d));
    }

    inline insert_return_type insert(const std::pair<key_type, mapped_type> t)
    {
        auto hash = hasher(t.first);
        size_t ttl = Ext::tab(hash, 0);

        auto result = Base_t::insert(t);
        if (result.second)
        {
            auto currsize = ++ll_elem[ttl];
            if (currsize > ll_thresh[ttl]) growTab(ttl);
        }
        return result;
    }

    size_t erase(const key_type k)
    {
        auto hash     = hasher(k);
        size_t ttl    = Ext::tab(hash, 0);
        size_t nk     = Base_t::erase(k);
        ll_elem[ttl] -= nk;
        return nk;
    }



    std::pair<size_t, Bucket_t*> getTable(size_t i)
    {
        return (i < tl) ? std::make_pair(ll_size[i], ll_tab[i].get())
                        : std::make_pair(0,nullptr);
    }

    iterator begin()
    {
        auto temp = make_iterator(&ll_tab[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

    const_iterator cbegin() const
    {
        auto temp = make_citerator(&ll_tab[0][0].elements[0]);
        if (! temp->first) temp++;
        return temp;
    }

private:
    // Functions for finding buckets *******************************************

    inline void getBuckets(Hashed_t h, Bucket_t** mem) const
    {
        for (size_t i = 0; i < nh; ++i)
            mem[i] = getBucket(h, i);
    }

    inline Bucket_t* getBucket (Hashed_t h, size_t i) const
    {
        size_t tab = Ext::tab(h,0);
        return &(ll_tab[tab][Ext::loc(h,i)*ll_factor[tab]]);
    }



    // Size changes (GROWING) **************************************************
    void grow()
    {
        // necessary for base class not implemented since growing is triggered
        // individually for each subtable in this specialization
    }

    inline void growTab(size_t tab)
    {
        //TODO make this prettier
        if (grow_buffer.size()) return;
        size_t nsize   = std::floor(double(ll_elem[tab]) * alpha / double(bs));
        nsize          = std::max(nsize, ll_size[tab]+1);
        capacity      += nsize - ll_size[tab];
        double nfactor = double(nsize)      / fac_div;
        size_t nthresh = ll_elem[tab] * beta;

        auto ntable = std::make_unique<Bucket_t[]>(nsize);
        migrate(tab, ntable, nfactor);

        ll_tab[tab]    = std::move(ntable);
        ll_size[tab]   = nsize;
        ll_factor[tab] = nfactor;
        ll_thresh[tab] = nthresh;
        if (grow_buffer.size()) finalize_grow();

    }

    inline void migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, double nfactor)
    {
        size_t csize   = ll_size[tab];
        double cfactor = ll_factor[tab];
        for (size_t i = 0; i < csize; ++i)
        {
            Bucket_t& curr = ll_tab[tab][i];

            for (size_t j = 0; j < bs; ++j)
            {
                auto e    = curr.elements[j];
                if (! e.first) break;
                auto hash = hasher(e.first);

                for (size_t ti = 0; ti < nh; ++ti)
                {
                    if (i == size_t(Ext::loc(hash, ti) * cfactor))
                    {
                        if (! target[Ext::loc(hash, ti) * nfactor].insert(e))
                            grow_buffer.push_back(e);
                        break;
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
            Base_t::insert(e);
        }
        n = temp;
        std::vector<std::pair<key_type, mapped_type> > ttemp;
        std::swap(ttemp, grow_buffer);
    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF, class Conf>
class CuckooTraits<CuckooIndependent2L<K,D,HF,Conf> >
{
public:
    using Specialized_t  = CuckooIndependent2L<K,D,HF,Conf>;
    using Base_t         = CuckooMultiBase<Specialized_t>;
    using Config_t       = Conf;

    using key_type       = K;
    using mapped_type    = D;

    static constexpr size_t bs = Conf::bs;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t nh = Conf::nh;

    using Hasher_t       = Hasher<K, HF, ct_log(tl), nh, true, true>;
    using Bucket_t       = Bucket<K,D,bs>;
};



// Iterator increment **********************************************************

template<class K, class D, class HF, class Conf>
class iterator_incr<CuckooIndependent2L<K,D,HF,Conf> >
{
private:
    using ipointer = std::pair<const K,D>*;
    static constexpr size_t tl = Conf::tl;
    static constexpr size_t bs = Conf::bs;

public:
    using Table_t   = CuckooIndependent2L<K,D,HF,Conf>;

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
        size_t size = table.ll_size[tab] - 1;
        end_tab = reinterpret_cast<ipointer>(&table.ll_tab[tab][size].elements[bs-1]);
        return reinterpret_cast<ipointer>(&table.ll_tab[tab][0].elements[0]);
    }

    void initialize_tab(ipointer ptr)
    {
        for (size_t i = 0; i < tl; ++i)
        {
            ipointer tab_b_ptr = reinterpret_cast<ipointer>
                (&table.ll_tab[i][0].elements[0]);
            ipointer tab_e_ptr = tab_b_ptr + table.ll_size[i]*bs;

            if (tab_b_ptr <= ptr && ptr < tab_e_ptr)
            {
                tab = i;
                end_tab = tab_e_ptr;
                return;
            }
        }
    }

};
