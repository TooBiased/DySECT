#pragma once

#include "cuckoo_base.h"

template<class K0, class D0, class HF0, class Conf0>
class CuckooIndependentBase;



// *****************************************************************************
// IN PLACE GROWING ************************************************************
// *****************************************************************************

template<class K, class D, class HF = std::hash<K>,
         class Conf = Config<> >
class CuckooDeAmortized : public CuckooTraits<CuckooDeAmortized<K,D,HF,Conf> >::Base_t
{
private:
    using This_t         = CuckooDeAmortized<K,D,HF,Conf>;
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
    static constexpr size_type min_buckets  = 4096;
    static constexpr size_type grow_step    = 32;

public:
    CuckooDeAmortized(size_type cap = 0      , double size_constraint = 1.1,
                 size_type dis_steps = 0, size_type seed = 0)
        : Base_t(size_constraint, dis_steps, seed)
    {
        auto temp = reinterpret_cast<value_intern*>(operator new (max_size));
        table = std::unique_ptr<value_intern[]>(temp);

        // SET CAPACITY, BITMASKS, THRESHOLD
        size_type tcap = size_type(double(cap) * size_constraint / double(bs));
        size_type tbit = min_buckets << 1; while (tbit <= tcap) tbit <<= 1;
        bitmask_large = tbit - 1;
        bitmask_small = bitmask_large >> 1;

        size_type doub = tcap -bitmask_small-1;
        doub /= grow_step;
        doub *= grow_step;

        capacity      = (doub+bitmask_small+1)*bs;
        bucket_cutoff = doub+bitmask_small+1;
        grow_thresh   = size_type(double(capacity+grow_step*bs)/alpha);

        std::fill(table.get(), table.get()+capacity, value_intern());
    }

    CuckooDeAmortized(const CuckooDeAmortized&) = delete;
    CuckooDeAmortized& operator=(const CuckooDeAmortized&) = delete;

    CuckooDeAmortized(CuckooDeAmortized&&) = default;
    CuckooDeAmortized& operator=(CuckooDeAmortized&&) = default;

private:
    using Base_t::n;
    using Base_t::capacity;
    using Base_t::grow_thresh;
    using Base_t::alpha;
    using Base_t::hasher;

    size_type bucket_cutoff;
    size_type bitmask_large;
    size_type bitmask_small;

    std::unique_ptr<value_intern[]> table;

    using Base_t::make_iterator;
    using Base_t::make_citerator;

public:
    using Base_t::insert;

    std::pair<size_type, Bucket_t*> getTable(size_type i)
    {
        return (! i) ? std::make_pair(bucket_cutoff, table.get())
                     : std::make_pair(0,nullptr);
    }

    inline iterator begin()
    {
        auto temp = make_iterator(&table[0]);
        if (! temp->first) temp++;
        return temp;
    }

    inline const_iterator cbegin() const
    {
        auto temp = make_citerator(&table[0]);
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
        size_type l = Ext::loc(h,i) & bitmask_large;
        if (l >= bucket_cutoff) l &= bitmask_small;
        return reinterpret_cast<Bucket_t*>(&(table[l*bs]));
    }



    // Size changes (GROWING) **************************************************

    inline void grow()
    {
        // std::cout << "cap:" << capacity
        //           << " cut:" << bucket_cutoff << std::endl;
        size_type ncap    = capacity + grow_step*bs;
        size_type ncutoff = bucket_cutoff + grow_step;
        grow_thresh       = size_type(double(ncap+grow_step*bs)/alpha);

        std::fill(table.get()+capacity, table.get()+ncap, value_intern());

        migrate(capacity, ncap);

        //auto ocap = capacity;
        //auto ocut = bucket_cutoff;
        capacity      = ncap;
        bucket_cutoff = ncutoff;

        //for (size_t i = gable.get)

        // size_type bla = 0;
        // size_type blu = 0;
        // for (size_type i = 0; i < capacity; ++i)
        // {
        //     auto curr = table[i];
        //     if (!curr.first) continue;
        //     ++bla;
        //     if (Base_t::find(curr.first) == Base_t::end())
        //     {
        //         //std::cout << i << std::endl;
        //         ++blu;
        //     }
        // }
        // std::cout << "n " << n << "  moved " << bla << "  lost " << blu << std::endl;

        if (ncutoff >= bitmask_large)
        {
            bitmask_small = bitmask_large;
            bitmask_large = (bitmask_large<<1) + 1;
            // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }

    }

    inline void migrate(size_type , size_type)
    {
        size_type flag  = bitmask_large  & (~bitmask_small);
        size_type ptr0  = (bucket_cutoff & bitmask_small) * bs;

        //static size_type counter = 0;

        //for (size_type i = 0; i < grow_step*bs; i+=bs)
        for (size_type i = 0; i < grow_step; ++i)
        {
            // Bucket_t* bucket0_ptr =
            //     reinterpret_cast<Bucket_t*>(table.get() + ((cap/bs) & bitmask_small) + i);
            // Bucket_t* bucket1_ptr =
            //     reinterpret_cast<Bucket_t*>(table.get() + capacity +i);
            //if (((i/bs)&bitmask_small)%bs) std::cout << "9222222222222222222222222222222222222222222222222222" << std::endl;
            Bucket_t* bucket0_ptr =
                reinterpret_cast<Bucket_t*>(table.get()+ptr0+i*bs);
            Bucket_t* bucket1_ptr =
                reinterpret_cast<Bucket_t*>(table.get()+capacity+i*bs);

            //static size_type fuck = 0;
            // if (fuck < 10)
            // {
            //     std::cout << capacity << " "
            //               << off0s+i << " "< bitmask_small << " "
            //               << bitmask_large << std::endl;
            //     fuck++;
            // }
            // for (size_type k = 0; k < bs; ++k)
            // {
            //     if (bucket1_ptr->elements[k].first) std::cout << "fuck " << std::flush;
            // }
            size_type c0 = 0;
            size_type c1 = 0;
            for (size_type j = 0; j < bs; ++j)
            {
                auto curr = bucket0_ptr->elements[j];
                bucket0_ptr->elements[j] = value_intern();
                if (!curr.first) break;

                Bucket_t* targets[nh];
                Hashed_t  hash = hasher(curr.first);
                getBuckets(hash, targets);

                // bool wtf = false;
                for (size_type t = 0; t < nh; ++t)
                {
                    if (targets[t] == bucket0_ptr)
                    {
                        if (Ext::loc(hash, t) & flag)
                        {
                            bucket1_ptr->elements[c1++] = curr;
                            // if (bucket1_ptr != reinterpret_cast<Bucket_t*>(table.get() + (Ext::loc(hash, t) & bitmask_large)*bs))
                            //     std::cout << "i:" << (i-capacity)/bs << "k" << std::endl;
                        }
                        else
                        {
                            bucket0_ptr->elements[c0++] = curr;
                            // if (bucket0_ptr != reinterpret_cast<Bucket_t*>(table.get() + (Ext::loc(hash, t) & bitmask_large)*bs))
                            //     std::cout << "g" << std::endl;
                        }
                        //wtf = true;
                        break;
                    }
                }

                // if (!wtf)
                // {
                //     std::cout <<   "i:" << (i - capacity)/bs
                //               << " b0:" << bucket0_ptr
                //               << " b1:" << bucket1_ptr
                //               << " stuff:" << std::flush;
                //     for (size_t z = 0; z < nh; ++z)
                //     {
                //         std::cout << " " << targets[z];
                //     }
                //     std::cout << std::endl;
                //     counter++;
                // }
            }
        }
        // if (counter)
        // { std::cout << counter << std::endl; }
    }
};



// Traits class defining types *************************************************

template<class K, class D, class HF,
         class Conf>
class CuckooTraits<CuckooDeAmortized<K,D,HF,Conf> >
{
public:
    using Specialized_t = CuckooDeAmortized<K,D,HF,Conf>;
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
class iterator_incr<CuckooDeAmortized<K,D,HF,Conf> >
{
public:
    using Table_t   = CuckooDeAmortized<K,D,HF,Conf>;
private:
    using size_type = typename Table_t::size_type;
    using pointer   = std::pair<const K,D>*;
    static constexpr size_type bs = Conf::bs;

public:
    iterator_incr(const Table_t& table_)
        : end_ptr(reinterpret_cast<pointer>(&table_.table[table_.capacity -1]))
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
