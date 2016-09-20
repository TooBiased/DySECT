#pragma once


#include <functional>
#include <cmath>
#include <memory>
#include <iostream>
#include <vector>
#include <tuple>
#include <random>

#define MIN_LLS 256

template <class T>
inline void print(std::ostream& out, const T& t, size_t w)
{
    out.width(w);
    out << t << " " << std::flush;
}



template <class Parent>
class dstrat_triv
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using HashSplitter = typename Parent::HashSplitter;

    dstrat_triv(Parent&, size_t, size_t) {}

    size_t  steps = 0;
    size_t* hist  = nullptr;
    
    bool insert(std::pair<Key, Data>, HashSplitter)
    {   return false;   }
};



template<class, class, size_t>
class Bucket;

template<class K, class D, class H = std::hash<K>,
         template<class> class DS = dstrat_triv,
         size_t TL = 128, size_t BS = 4>   // ALPHA (SIZE CONSTRAINT COULD BE PARAMETER BUT INTEGRAL TYPE PROBLEM)
class SpaceGrow
{
public:
    using Key       = K;
    using Data      = D;
    using FRet      = std::pair<bool, Data>;    

    SpaceGrow(size_t _capacity = 0, double size_constraint = 1.1,
              size_t dis_steps = 0, size_t seed = 0);

    ~SpaceGrow() = default;

    SpaceGrow(const SpaceGrow&) = delete;
    SpaceGrow& operator=(const SpaceGrow&) = delete;

    SpaceGrow(SpaceGrow&& rhs) = default;
    SpaceGrow& operator=(SpaceGrow&& rhs) = default;

    bool insert(Key k, Data d);
    bool insert(std::pair<Key, Data> t);
    FRet find  (Key k);
    bool remove(Key k);

    void printDist(std::ostream& out);
    void printHist(std::ostream& out);
    
private:    
    using  This_t     = SpaceGrow<K,D,H,DS,TL,BS>;
    using  Bucket_t   = Bucket<K,D,BS>;
    using  HashFct_t  = H;
    using  DisStrat_t = DS<This_t>;
    friend DisStrat_t;

    static constexpr size_t log(size_t k)
    {
        return (k) ? 1+log(k>>1) : 0;
    }
    
    union HashSplitter {
        std::uint64_t hash;
        struct
        {
            uint64_t tab1 : log(TL)-1;
            uint64_t tab2 : log(TL)-1;
            uint64_t loc1 : 32-log(TL)+1;
            uint64_t loc2 : 32-log(TL)+1;
        };
    };
    
    static_assert( sizeof(HashSplitter)==8,
                   "HashSplitter must be 64bit!" );
    
    //static_assert( !(((1<<log(TL))-1) & TL),
    //               "TL must be a power of two >0!");
    
    HashSplitter h(Key k)
    {
        HashSplitter a;
        a.hash = hasher(k);
        return a;
    }
    
public: //temporary should be removed

    size_t       nElements;
    size_t       curGrowAmount;
    size_t       curGrowTable;
    size_t       capacity;
    double       alpha;
    const size_t bs = BS;
    HashFct_t    hasher;
    DisStrat_t   displacer;
    
    alignas(64) size_t                      llb[TL];
    alignas(64) std::unique_ptr<Bucket_t[]> llt[TL];  // lower level tables

    void incElements();
    void grow();
    void migrate(size_t i, std::unique_ptr<Bucket_t[]>& target, size_t tBitmask);

    Bucket_t* getBucket1(HashSplitter h)
    {   return &(llt[h.tab1][(h.loc1 & llb[h.tab1])]);   }
    Bucket_t* getBucket2(HashSplitter h)
    {   return &(llt[h.tab2][(h.loc2 & llb[h.tab2])]);   }
    
};


template<class K, class D, size_t BS = 4>
class Bucket
{
public:
    using Key = K;
    using Data = D;
    using FRet = std::pair<bool, Data>;

    Bucket() { for (size_t i = 0; i < BS; ++i) elements[i] = std::make_pair(Key(), Data());}
    Bucket(const Bucket& rhs) = default;
    Bucket& operator=(const Bucket& rhs) = default;
    
    bool   insert(Key k, Data d);
    bool   insert(std::pair<Key, Data> t);
    FRet   find  (Key k);
    bool   remove(Key k);
    FRet   pop   (Key k);
    
    int    probe (Key k);

    bool   space ();
    std::pair<Key, Data> get(size_t i);
    std::pair<Key, Data> replace(size_t i, std::pair<Key, Data> t);
    
    std::pair<Key, Data> elements[BS];
};




template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
SpaceGrow<K,D,HF,DS,TL,BS>::SpaceGrow(size_t _capacity, double size_constraint,
                                      size_t dis_steps, size_t seed)
    : nElements(0), alpha(size_constraint), displacer(*this, dis_steps, seed)
{
    std::cout << "this_size:      " << sizeof(This_t)       << std::endl;
    std::cout << "dis_strat_size: " << sizeof(DisStrat_t)   << std::endl;
    std::cout << "hasher_size:    " << sizeof(HashFct_t)    << std::endl;
        
    double dIni = double(_capacity) * size_constraint / double(TL);
    if (dIni < MIN_LLS)
    {
        for (size_t i = 0; i < TL; ++i)
        {
            llb[i] = (MIN_LLS/BS) - 1;
            llt[i] = std::make_unique<Bucket_t[]>(MIN_LLS/BS);
        }
        curGrowAmount = MIN_LLS;
        curGrowTable  = 0;
        capacity      = MIN_LLS*TL;
            
    }
    else
    {
        size_t iIni = MIN_LLS;
        while (dIni > (iIni << 1)) iIni <<= 1;

        size_t gIni = std::floor(double(_capacity) * size_constraint / double(iIni))-TL;;
            
        for (size_t i = gIni; i < TL; ++i)
        {
            llb[i] = (iIni/BS)-1;
            llt[i] = std::make_unique<Bucket_t[]>(iIni/BS);
        }

        curGrowAmount = iIni;
        curGrowTable  = gIni;
        capacity       = (gIni+TL) * iIni;
            
        iIni         <<= 1;

        for (size_t i = 0; i < gIni; ++i)
        {
            llb[i] = (iIni/BS)-1;
            llt[i] = std::make_unique<Bucket_t[]>(iIni/BS);
        }

    }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(Key k, Data d)
{
    return insert(std::make_pair(k,d));
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(std::pair<Key, Data> t)
{
    auto hash = h(t.first);
        
    auto p1 = getBucket1(hash)->probe(t.first);//llt[ hash.tab1 ].probe( k, hash.loc1 );
    auto p2 = getBucket2(hash)->probe(t.first);//llt[ hash.tab2 ].probe( k, hash.loc2 );
    
    if ((p1 < 0) || (p2 < 0)) return false;

    auto r = false;
    
    // SHOULD CHECK, IF ALREADY INCLUDED
    if (p1 > p2)
    {
        // insert into p1
        r = getBucket1(hash)->insert(t);
    }
    else if (p2 > p1)
    {
        // insert into p2
        r = getBucket2(hash)->insert(t);
    }
    else if (p1)
    {
        // insert into p1
        r = getBucket1(hash)->insert(t);
    }
    else
    {
        // no space => displace stuff
        r = displacer.insert(t, hash);
    }
    
    if (r) incElements();

    return r;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
typename SpaceGrow<K,D,HF,DS,TL,BS>::FRet SpaceGrow<K,D,HF,DS,TL,BS>::find(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k);

    if (p1.first) return p1;
    else          return getBucket2(hash)->find(k);
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k);
    
    if (p1) return true;
    else    return getBucket2(hash)->remove(k);
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::grow()
{
    size_t nsize = (curGrowAmount << 1)/BS; //double the previous size
    auto   ntab  = std::make_unique<Bucket_t[]>(nsize);
    migrate(curGrowTable, ntab, nsize-1);
    llb[curGrowTable] = nsize-1;
    llt[curGrowTable] = std::move(ntab);
    capacity          += curGrowAmount;
    if (++curGrowTable == TL) { curGrowTable = 0; curGrowAmount <<= 1; }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::incElements()
{
    ++nElements;
    if (capacity + curGrowAmount < nElements * alpha) grow();
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::migrate(size_t tab, std::unique_ptr<Bucket_t[]>& target, size_t bitmask)
{
    for (size_t i = 0; i <= llb[tab]; ++i) //Bucket_t* curr = &(llt[i][0]); curr <= &(llt[i][llb[i]]); ++curr)
    {
        Bucket_t* curr = &(llt[tab][i]);
        for (size_t j = 0; j < BS; ++j)
        {
            auto e    = curr->elements[j];
            if (! e.first) break;
            auto hash = h(e.first);
            if      (getBucket1(hash) == curr) target[bitmask & hash.loc1].insert(e.first, e.second);
            else if (getBucket2(hash) == curr) target[bitmask & hash.loc2].insert(e.first, e.second);
            else { std::cout << "something is wrong neither in first, nor second bucket." << std::endl; }
        }
    }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::printDist(std::ostream& out)
{
    size_t gHist[BS];
    for (size_t i = 0; i < BS; ++i) gHist[i] = 0;

    print (out, "# tab", 5);
    print (out, "full" , 6);
    print (out, "1 sp" , 6);
    print (out, "2 sp" , 6);
    print (out, "empt" , 6);
    print (out, "n"    , 8);
    print (out, "cap"  , 8);
    out << std::endl;
    
    for (size_t tl = 0; tl < TL; ++tl)
    {
        size_t lHist[BS];
        for (size_t i = 0; i < BS; ++i) lHist[i] = 0;
        
        for (size_t j = 0; j <= llb[tl]; ++j)
        {
            auto a = llt[tl][j].probe(0);
            if (a >= 0) ++lHist[a];
        }

        size_t n = 0;
        print (out, tl, 5);
        for (size_t i = 0; i < BS; ++i)
        {
            print (out, lHist[i], 6);
            n += lHist[i] * (BS - i);
            gHist[i] += lHist[i];
        }
        print (out, n, 8);
        print (out, (llb[tl]+1)*BS, 8);
        out << std::endl;
    }

    size_t n = 0;
    print (out, "#all", 5);
    for (size_t i = 0; i < BS; ++i)
    {
        print (out, gHist[i], 6);
        n += gHist[i] * (BS - i);
    }
    print (out, n, 8);
    print (out, curGrowAmount * (TL+curGrowTable), 8);
    out << std::endl;    
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::printHist(std::ostream& out)
{
    print(out, "# steps", 7);
    print(out, "nFitted", 8);
    out << std::endl;
    
    for (size_t i = 0; i < displacer.steps; ++i)
    {
        print(out, i, 7);
        print(out, displacer.hist[i], 8);
        out << std::endl;
    }
}




template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::insert(Key k, Data d)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first) continue;
        elements[i].first  = k;
        elements[i].second = d;

        return true;
    }

    return false;
}

template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::insert(std::pair<Key,Data> t)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first) continue;
        
        elements[i]  = t;
        return true;
    }

    return false;
}

template<class K, class D, size_t BS>
typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::find(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first )      return std::make_pair(false, Data());
        if ( elements[i].first == k ) return std::make_pair(true, elements[i].second);
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::remove(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            for (size_t j = i+1; j < BS; ++j)
            {
                if (elements[j].first) { elements[i] = elements[j]; i = j; }
                else break;
            }
            elements[i] = std::make_pair(Key(), Data());
            return true;                                    
        }
        else if (! elements[i].first)
        {
            break;
        }
    }
    return false;
}

template<class K, class D, size_t BS>
typename Bucket<K,D,BS>::FRet Bucket<K,D,BS>::pop(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (elements[i].first == k)
        {
            Data d = elements[i].second;
            for (size_t j = i+1; j < BS; ++j)
            {
                if (elements[j].first) { elements[i] = elements[j]; i = j; }
                else break;
            }
            elements[i] = std::make_pair(Key(), Data());
            return std::make_pair(true, d);                                    
        }
        else if (! elements[i].first)
        {
            break;
        }
    }
    return std::make_pair(false, Data());
}

template<class K, class D, size_t BS>
int Bucket<K,D,BS>::probe(Key k)
{
    for (size_t i = 0; i < BS; ++i)
    {
        if (!elements[i].first)      return BS - i;
        if ( elements[i].first == k) return -1;
        
    }
    return 0;
}

template<class K, class D, size_t BS>
bool Bucket<K,D,BS>::space()
{
    return !elements[BS-1].first;
}

template<class K, class D, size_t BS>
std::pair<K, D> Bucket<K,D,BS>::get(size_t i)
{
    return elements[i];
}

template<class K, class D, size_t BS>
std::pair<K, D> Bucket<K,D,BS>::replace(size_t i, std::pair<K, D> newE)
{
    auto temp = elements[i];
    elements[i] = newE;
    return temp;
}
