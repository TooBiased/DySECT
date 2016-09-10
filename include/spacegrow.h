#pragma once


#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

#define MIN_LLS 256



template<class Parent>
class dstrat_bfs
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using Parent_t     = typename Parent::This_t;
    using HashSplitter = typename Parent::HashSplitter;
    using LLTable_t    = typename Parent::LLTable_t;
    using Bucket_t     = typename Parent::Bucket_t;

    using BFSQueue     = std::vector<std::tuple<Key, size_t, size_t> >;

    static constexpr size_t max_vis_nodes = 1024;
    
    inline static bool expand(Parent_t& p, BFSQueue& q, size_t tab, size_t loc)
    {
        Bucket_t tbuck = p.llt[tab].getBucket(loc);
        if (!tbuck[4-1].first) return true;

        for (size_t i = 0; i < 4; ++i)
        {
            Key k = tbuck[i].first;
            auto hash = p.h(k);
            if ((tab == hash.tab1) && (loc == hash.loc1))
            {
                q.push_back(std::tuple<Key, size_t, size_t>(k, hash.tab2, hash.loc2));
            } else
            {
                q.push_back(std::tuple<Key, size_t, size_t>(k, hash.tab1, hash.loc1));
            }
        }
        return false;
    }

    inline static bool rollBackDisplacements(Parent_t& p, BFSQueue& q, size_t i)
    {
        return false;
    }
    
    inline static bool insert(Parent_t& p, Key k, Data , HashSplitter hash)
    {
        std::vector<std::tuple<Key, size_t, size_t> > bq;
        bq.push_back(std::tuple<Key, size_t, size_t>(k, hash.tab1, hash.loc1));
        bq.push_back(std::tuple<Key, size_t, size_t>(k, hash.tab2, hash.loc2));
        for (size_t i = 0; i < max_vis_nodes; ++i)
        {
            Key k;
            size_t tab;
            size_t loc;
            std::tie(k, tab, loc) = bq[i];
            if (expand(p, bq, tab, loc)) return rollBackDisplacements(p, bq, i);
        }
        
        return false;
    }
};


template<class, class, size_t>
class Bucket;

template<class, class, size_t>
class LowerLevelTable;

template<class K, class D, class H = std::hash<K>,
         template<class> class DS = dstrat_bfs,
         size_t TL = 256, size_t BS = 4>   // ALPHA (SIZE CONSTRAINT COULD BE PARAMETER BUT INTEGRAL TYPE PROBLEM)
class SpaceGrow
{
public:
    using Key       = K;
    using Data      = D;
    using FRet      = std::pair<bool, Data>;
    using HashFct   = H;
    

    SpaceGrow(size_t capacity, double size_constraint = 1.05)
    {
        double dIni = double(capacity) * size_constraint / double(TL);
        if (dIni < MIN_LLS)
        {
            for (size_t i = 0; i < TL; ++i)
            {
                //std::cout << i <<  ":" << std::flush;
                lls[i] = MIN_LLS;
                llt[i] = LLTable_t(MIN_LLS);
            }
        } else
        {
            size_t iIni = MIN_LLS;
            while (dIni > (iIni << 1)) iIni <<= 1;

            size_t gIni = std::ceil(double(capacity) * size_constraint / double(iIni))-TL;;
            
            for (size_t i = gIni; i < TL; ++i)
            {
                //std::cout << i << ":" << std::flush;
                lls[i] = iIni;
                llt[i] = LLTable_t(iIni);
            }

            iIni <<= 1;

            for (size_t i = 0; i < gIni; ++i)
            {
                //std::cout << i <<  ":" << std::flush;
                lls[i] = iIni;
                llt[i] = LLTable_t(iIni);
            }

            std::cout << "space: " << (gIni+TL)*(iIni >> 1) << "cells" << std::endl;
        }
    }

    ~SpaceGrow() = default;

    SpaceGrow(const SpaceGrow&) = delete;
    SpaceGrow& operator=(const SpaceGrow&) = delete;

    SpaceGrow(SpaceGrow&& rhs) = default;
    SpaceGrow& operator=(SpaceGrow&& rhs) = default;

    bool insert(Key k, Data d);
    FRet find  (Key k);
    bool remove(Key k);
    
private:
    using  This_t     = SpaceGrow<K,D,H,DS,TL,BS>;
    using  LLTable_t  = LowerLevelTable<K,D,BS>;
    using  Bucket_t   = Bucket<K,D,BS>;
    using  DisStrat_t = DS<This_t>;
    friend DisStrat_t;
    //constexpr size_t nll = TL; // number of lower level tables
    
    size_t    lls[TL];  // sizes of lower level tables
    LLTable_t llt[TL];  // lower level tables

    HashFct   hasher;

    union HashSplitter {
        std::uint64_t hash;
        struct
        {
            uint64_t tab1: 8;
            uint64_t tab2: 8;
            uint64_t loc1: 24;
            uint64_t loc2: 24;
        };
    };
    static_assert( sizeof(HashSplitter)==8,
                   "Size of HashSep Object does is not 64bit=8Byte!" );

    HashSplitter h(Key k)
    {
        HashSplitter a;
        a.hash = hasher(k);
        return a;
    }
};



template<class K, class D, size_t BS = 4>
class LowerLevelTable
{
public:
    using Key  = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;
    using Bucket_t = Bucket<K, D, BS>;

    LowerLevelTable() : btable(nullptr) { } // has to have default cons will fail if table is used
    
    LowerLevelTable(size_t actual_size_not_capacity)
          // Has to be initialized with 2**k greater than MIN_LLS
        : bitmask((actual_size_not_capacity/BS)-1)
    {
        btable = new Bucket_t[bitmask+1];
        for (size_t i = 0; i < bitmask+1; ++i) { btable[i] = Bucket_t(); }
    }
    ~LowerLevelTable()
    {   delete[] btable;   }

    LowerLevelTable(const LowerLevelTable&) = delete;
    LowerLevelTable& operator=(const LowerLevelTable&) = delete;

    LowerLevelTable(LowerLevelTable&& rhs) : bitmask(rhs.bitmask), btable(nullptr)
    {   std::swap(btable, rhs.btable);   }
    LowerLevelTable& operator=(LowerLevelTable&& rhs)
    {   bitmask = rhs.bitmask; std::swap(btable, rhs.btable); return *this;  }

    bool   insert(Key k, Data d, size_t loc);
    FRet   find  (Key k,         size_t loc);
    bool   remove(Key k,         size_t loc);

    int    probe (Key k,         size_t loc);
    Bucket_t getBucket(size_t loc) { return btable[bucket(loc)]; }
    
private:
    //constexpr size_t sB = BS;

    size_t bitmask;
    Bucket_t* btable;

    size_t bucket(size_t loc) { return loc & bitmask; };
};



template<class K, class D, size_t BS = 4>
class Bucket
{
public:
    using Key = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    Bucket() { for (size_t i = 0; i < BS; ++i) p[i] = std::make_pair(Key(), Data());}
    Bucket(const Bucket& rhs) = default;
    Bucket& operator=(const Bucket& rhs) = default;
    
    bool   insert(Key k, Data d);
    FRet   find  (Key k);
    bool   remove(Key k);
    
    int    probe (Key k);
    
    
private:
    std::pair<Key, Data> p[BS];
};






template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(Key k, Data d)
{
    auto hash = h(k);
        
    auto p1 = llt[ hash.tab1 ].probe( k, hash.loc1 );
    auto p2 = llt[ hash.tab2 ].probe( k, hash.loc2 );

    if ((p1 < 0) || (p2 < 0)) return false;
    
    // SHOULD CHECK, IF ALREADY INCLUDED
    if (p1 > p2)
    {
        // insert into p1
        auto b = llt[ hash.tab1 ].insert(k, d, hash.loc1);
        return b;
    }
    else if (p2 > p1)
    {
        // insert into p2
        auto b = llt[ hash.tab2 ].insert(k, d, hash.loc2);
        return b;
    }
    else if (p1)
    {
        // insert into p1
        auto b = llt[ hash.tab1 ].insert(k, d, hash.loc1);
        return b;
    }
    else
    {
        DisStrat_t::insert(*this, k, d, hash);
        // no space, make space
        //std::cout << "no space left, make space" << std::endl;
        return false;
    }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
typename SpaceGrow<K,D,HF,DS,TL,BS>::FRet SpaceGrow<K,D,HF,DS,TL,BS>::find(Key k)
{
    auto hash = h(k);
    auto p1 = llt[hash.tab1].find(k, hash.loc1);

    if (p1.first) return p1;
    else          return llt[hash.tab2].find(k, hash.loc2);
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = llt[hash.tab1].remove(k, hash.loc1);
    if (p1) return true;
    else    return llt[hash.tab2].remove(k,hash.loc2);
}







template<class K, class D, size_t BS>
bool LowerLevelTable<K,D,BS>::insert(Key k, Data d, size_t loc)
{
    return btable[bucket(loc)].insert(k,d);
}

template<class K, class D, size_t BS>
typename LowerLevelTable<K,D,BS>::FRet LowerLevelTable<K,D,BS>::find(Key k, size_t loc)
{
    return btable[bucket(loc)].find(k);
}

template<class K, class D, size_t BS>
bool LowerLevelTable<K,D,BS>::remove(Key k, size_t loc)
{
    return btable[bucket(loc)].remove(k);
}

template<class K, class D, size_t BS>
int LowerLevelTable<K,D,BS>::probe(Key k, size_t loc)
{
    //std::cout << "loc:" << loc << "  bucket:" << bucket(loc) << std::endl;
    return btable[bucket(loc)].probe(k);
}





template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::insert(Key k, Data d)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (p[i].first) continue;
        p[i].first  = k;
        p[i].second = d;

        return true;
    }

    return false;
}

template<class K, class D, size_t BS>
typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::find(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!p[i].first )      return std::make_pair(false, Data());
        if ( p[i].first == k ) return std::make_pair(true, p[i].data);
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::remove(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (p[i].first == k)
        {
            for (size_t j = i+1; j < BS; ++j)
            {
                if (p[j].first) { p[i] = p[j]; i = j; }
                else break;
            }
            p[i] = std::make_pair(Key(), Data());
            return true;                                    
        }
        else if (! p[i].first)
        {
            return false;
        }
    }
}

template<class K, class D, size_t BS>
int Bucket<K,D,BS>::probe(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!p[i].first)      return BS - i;
        if ( p[i].first == k) return -1;
        
    }
    return 0;
}
