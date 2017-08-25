#pragma once

/*******************************************************************************
 * include/cuckoo_base.h
 *
 * CuckooMultiBase implements the basics of d-ary bucket cuckoo
 * hashing.  Implemented functions include insert, find, erase ... .
 * Inheriting classes only have to implement the get bucket functions
 * (and Specialize CuckooTraits, iterator_incrr).  CRTP is used to
 * eliminate vtable lookups.
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
#include <limits>

#include "bucket.h"
#include "hasher.h"
#include "iterator_base.h"
#include "displacement_strategies/main_strategies.h"

// CRTP base class for all cuckoo tables, this encapsulates
// main cuckoo table functionality (insert, find, and remove)

namespace dysect
{

    template<class T>
    class cuckoo_traits;
/* EXAMPLE IMPLEMENTATION
    {
    public:
        using Specialized_t  = T;
        using Base_t         = cuckoo_base<T>;
        using Config_t       = cuckoo_config<...>;

        using key_type       = ... ;
        using mapped_type    = ... ;

        static constexpr size_t tl = ... ;
        static constexpr size_t bs = ... ;
        static constexpr size_t nh = ... ;

        union hasher_type    = hasher<key_type, HashFct, ...>;
        using bucket_type    = bucket<key_type, mapped_type, bs>;
    };*/

    template<class T>
    class iterator_incr;

    class hist_count
    {
    public:
        hist_count(size_t s) : steps(s), hist(new size_t[s])
        { for (size_t i = 0; i < s; ++i) { hist[i] = 0; } }

        void add(size_t i) { auto ind = (i<steps) ? i:steps-1; ++hist[ind];}

        const size_t steps;
        std::unique_ptr<size_t[]> hist;
    };

    class no_hist_count
    {
    public:
        no_hist_count(size_t = 0) { }
        void add(size_t) { }
        static constexpr size_t  steps = 0;
        static constexpr size_t* hist  = nullptr;
    };

    template<size_t BS = 8, size_t NH = 3, size_t TL = 256,
             template <class> class DisStrat = cuckoo_displacement::trivial,
             class HistCount = no_hist_count>
    struct cuckoo_config
    {
        static constexpr size_t bs = BS;
        static constexpr size_t tl = TL;
        static constexpr size_t nh = NH;
        static constexpr size_t sbs = 4;

        template <class T>
        using dis_strat_type  = DisStrat<T>;

        using hist_count_type = HistCount;
    };





    template<class SCuckoo>
    class cuckoo_base
    {
    private:
        using  this_type        = cuckoo_base<SCuckoo>;
        using  specialized_type = typename cuckoo_traits<SCuckoo>::specialized_type;
        using  dis_strat_type   = typename cuckoo_traits<SCuckoo>::config_type::template dis_strat_type<this_type>;
        using  hist_count_type  = typename cuckoo_traits<SCuckoo>::config_type::hist_count_type;
        using  bucket_type      = typename cuckoo_traits<SCuckoo>::bucket_type;
        using  hasher_type      = typename cuckoo_traits<SCuckoo>::hasher_type;
        using  hashed_type      = typename hasher_type::hashed_type;

        friend specialized_type;
        friend dis_strat_type;

    public:
        using key_type        = typename cuckoo_traits<SCuckoo>::key_type;
        using mapped_type     = typename cuckoo_traits<SCuckoo>::mapped_type;
        using value_type      = std::pair<const key_type, mapped_type>;
        using iterator        = iterator_base<iterator_incr<specialized_type> >;
        using const_iterator  = iterator_base<iterator_incr<specialized_type>, true>;
        using size_type       = size_t;
        using difference_type = std::ptrdiff_t;
        // using hasher          = Hash;
        // using key_equal       = KeyEqual
        // using allocator_type  = Allocator
        using reference       = value_type&;
        using const_reference = const value_type&;
        // using pointer         = std::allocator_traits<Allocator>::pointer;
        // using const_pointer   = std::allocator_traits<Allocator>::const_pointer;
        using insert_return_type   = std::pair<iterator, bool>;

        using local_iterator       = void;
        using const_local_iterator = void;
        using node_type            = void;
    private:
        using  value_intern    = std::pair<key_type, mapped_type>;

    public:
        cuckoo_base(double size_constraint = 1.1,
                    size_type dis_steps = 0, size_type seed = 0);
        ~cuckoo_base() = default;
        cuckoo_base(const cuckoo_base&     ) = delete;
        cuckoo_base(      cuckoo_base&& rhs);
        cuckoo_base& operator=(const cuckoo_base&     ) = delete;
        cuckoo_base& operator=(      cuckoo_base&& rhs)
            {
                n = rhs.n; capacity = rhs.capacity;
                grow_thresh = rhs.grow_thresh; alpha = rhs.alpha;
                return *this;
            }

    private:
    // Members *****************************************************************
        size_type       n;
        size_type       capacity;
        size_type       grow_thresh;
        double          alpha;
        hasher_type     hasher;
        dis_strat_type  displacer;
        hist_count_type hcounter;
        static constexpr size_type bs = cuckoo_traits<specialized_type>::bs;
        static constexpr size_type tl = cuckoo_traits<specialized_type>::tl;
        static constexpr size_type nh = cuckoo_traits<specialized_type>::nh;

    public:
    // Basic Hash Table Functionality ******************************************
        iterator              find  (const key_type& k);
        const_iterator        find  (const key_type& k) const;
        insert_return_type    insert(const key_type& k, const mapped_type& d);
        insert_return_type    insert(const value_intern& t);
        size_type             erase (const key_type& k);

    // Easy use Accessors for std compliance ***********************************
        inline iterator       begin ();       // see specialized_type
        inline const_iterator begin () const { return static_cast<specialized_type*>(this)->cbegin(); }
        inline const_iterator cbegin() const; // see specialized_type
        inline iterator       end   ()       { return make_iterator(nullptr); }
        inline const_iterator end   () const { return static_cast<specialized_type*>(this)->cend(); }
        inline const_iterator cend  () const { return make_citerator(nullptr); }

        mapped_type&          at    (const key_type& k);
        const mapped_type&    at    (const key_type& k) const;
        mapped_type&          operator[](const key_type& k);
        size_type             count (const key_type& k) const;

    // Global fill state *******************************************************
        inline size_type      empty()    const { return (n == 0); }
        inline size_type      size()     const { return n; }
        inline size_type      max_size() const { return (1ull << 32)*bs; }

        inline void           clear()
            { auto temp = specialized_type(0, alpha); (*this) = temp; }


    private:
    // Easy iterators **********************************************************
        inline iterator       make_iterator (      value_intern* pos) const
            { return iterator      (pos, *static_cast<const specialized_type*>(this)); }

        inline const_iterator make_citerator(const value_intern* pos) const
            { return const_iterator(pos, *static_cast<const specialized_type*>(this)); }

    // implementation specific functions (static polymorph) ********************
        inline void           inc_n() { ++n; }
        inline void           dec_n() { --n; }
        inline void           get_buckets(hashed_type h, bucket_type** mem) const
            { return static_cast<const specialized_type*>(this)->get_buckets(h, mem); }
        inline bucket_type*   get_bucket (hashed_type h, size_type i) const
            { return static_cast<const specialized_type*>(this)->get_bucket(h, i); }

    public:
    // auxiliary functions for testing *****************************************
        void                  clearHist();
        void                  print_init_data(std::ostream& out);
        static void           print_init_header(std::ostream& out)
        {
                out.width(6); out << "bsize";
                out.width(6); out << "ntabl";
                out.width(6); out << "nhash";
                out.width(9); out << "f_cap";
                out << std::flush;
        }

        void explicit_grow()
        {
                static_cast<specialized_type*>(this)->grow();
        }
    };



// Constructors and Appointments ***********************************************

    template<class SCuckoo>
    cuckoo_base<SCuckoo>::cuckoo_base(double size_constraint,
                                              size_type dis_steps, size_type seed)
        : n(0), capacity(0), grow_thresh(std::numeric_limits<size_type>::max()),
          alpha(size_constraint),
          displacer(*this, dis_steps, seed),
          hcounter(dis_steps)
    { }

    template<class SCuckoo>
    cuckoo_base<SCuckoo>::cuckoo_base(cuckoo_base&& rhs)
        : n(rhs.n), capacity(rhs.capacity), alpha(rhs.alpha),
          displacer(*this, std::move(rhs.displacer))
    { }



// Implementation of main functionality ****************************************

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::iterator
    cuckoo_base<SCuckoo>::find(const key_type& k)
    {
        auto hash = hasher(k);

        for (size_type i = 0; i < nh; ++i)
        {
            bucket_type* tb = get_bucket(hash, i);
            value_intern*   tp = tb->find_ptr(k);
            if (tp) return make_iterator(tp);
        }
        return end();
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::const_iterator
    cuckoo_base<SCuckoo>::find(const key_type& k) const
    {
        auto hash = hasher(k);

        for (size_type i = 0; i < nh; ++i)
        {
            bucket_type* tb = get_bucket(hash, i);
            value_intern*   tp = tb->find_ptr(k);
            if (tp) return make_citerator(tp);
        }
        return end();
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::insert_return_type
    cuckoo_base<SCuckoo>::insert(const key_type& k, const mapped_type& d)
    {
        return insert(std::make_pair(k,d));
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::insert_return_type
    cuckoo_base<SCuckoo>::insert(const value_intern& t)
    {
        if (n > grow_thresh) static_cast<specialized_type*>(this)->grow();
        auto hash = hasher(t.first);

        std::pair<int,value_intern*> max = std::make_pair(0, nullptr);
        for (size_type i = 0; i < nh; ++i)
        {
            auto temp = get_bucket(hash, i)->probe_ptr(t.first);

            if (temp.first < 0)
                return std::make_pair(make_iterator(temp.second), false);
            max = (max.first > temp.first) ? max : temp;
        }

        if (max.first > 0)
        {
            *max.second = t;
            hcounter.add(0);
            static_cast<specialized_type*>(this)->inc_n();
            return std::make_pair(make_iterator(max.second), true);
        }

        int  srch = -1;
        value_intern* pos  = nullptr;
        std::tie(srch, pos) = displacer.insert(t, hash);
        if (srch >=0)
        {
            hcounter.add(srch);
            static_cast<specialized_type*>(this)->inc_n();
            return std::make_pair(make_iterator(pos), true);
        }

        return std::make_pair(end(), false);
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::size_type
    cuckoo_base<SCuckoo>::erase(const key_type& k)
    {
        auto hash = hasher(k);
        for (size_type i = 0; i < nh; ++i)
        {
            bucket_type* tb = get_bucket(hash, i);
            if (tb->remove(k))
            {
                static_cast<specialized_type*>(this)->dec_n();
                return 1;
            }
        }
        return 0;
    }



// Accessor Implementations ****************************************************

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::mapped_type&
    cuckoo_base<SCuckoo>::at(const key_type& k)
    {
        auto a = static_cast<specialized_type*>(this)->find(k);
        if (a == end()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class SCuckoo>
    inline const typename cuckoo_base<SCuckoo>::mapped_type&
    cuckoo_base<SCuckoo>::at(const key_type& k) const
    {
        auto a = static_cast<const specialized_type*>(this)->find(k);
        if (a == cend()) throw std::out_of_range("cannot find key");
        else return (*a).second;
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::mapped_type&
    cuckoo_base<SCuckoo>::operator[](const key_type& k)
    {
        auto t = static_cast<specialized_type*>(this)->insert(k, mapped_type());
        return (*t.first).second;
    }

    template<class SCuckoo>
    inline typename cuckoo_base<SCuckoo>::size_type
    cuckoo_base<SCuckoo>::count(const key_type& k) const
    {
        return (static_cast<const specialized_type*>(this)->find(k) != cend()) ? 1 : 0;
    }



// Print Parameter Functions ***************************************************

/*template<class SCuckoo>
  inline static void cuckoo_base<SCuckoo>::print_init_header(std::ostream& out)
  {
  out.width(6); out << "bsize";
  out.width(6); out << "ntabl";
  out.width(6); out << "nhash";
  out.width(9); out << "f_cap";
  out << std::flush;
  }*/

    template<class SCuckoo>
    inline void cuckoo_base<SCuckoo>::print_init_data(std::ostream& out)
    {
        out.width(6); out << bs;
        out.width(6); out << tl;
        out.width(6); out << nh;
        out.width(9); out << capacity;
        out << std::flush;
    }

    template<class SCuckoo>
    inline void cuckoo_base<SCuckoo>::clearHist()
    {
        for (size_type i = 0; i < hcounter.steps; ++i) hcounter.hist[i] = 0;
    }

} // namespace dysect
