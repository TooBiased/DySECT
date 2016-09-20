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
inline void out (const T& t, size_t w)
{
    std::cout.width(w);
    std::cout << t << " " << std::flush;
}

template<class Parent>
class dstrat_bfs
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using Parent_t     = typename Parent::This_t;
    using HashSplitter = typename Parent::HashSplitter;
    using Bucket_t     = typename Parent::Bucket_t;

    using BFSQueue     = std::vector<std::tuple<Key, int, Bucket_t*> >;

    Parent&      tab;
    const size_t visit;
    size_t*      hist;

    dstrat_bfs(Parent& parent, size_t steps = 256, size_t t = 0)
        : tab(parent), visit((steps) ? steps : 256)
    {
        (void)t; /* parameter t is for symmetry with "rwalk" therefore unused*/
        hist = new size_t[visit];
        for (size_t i = 0; i < visit; ++i)
        {   hist[i] = 0;   }
    }
    ~dstrat_bfs() { delete[] hist; }
    
    inline bool expand(BFSQueue& q, size_t index)
    {
        Bucket_t* b = std::get<2>(q[index]);

        if (q.size() >= visit) return false;
        
        for (size_t i = 0; i < tab.bs; ++i)
        {
            Key k = b->p[i].first;
            
            auto hash = tab.h(k);
            Bucket_t* b1 = tab.getBucket1(hash);
            Bucket_t* b2 = tab.getBucket2(hash);
            
            if        (b == b1)
            {
                q.emplace_back(k, index, b2);
                if (! b2->p[tab.bs-1].first) return true;
            }
            else if (b == b2)
            {
                q.emplace_back(k, index, b1);
                if (! b1->p[tab.bs-1].first) return true;
            }
            else
            {
                std::cout << "unexpectedly in wrong bucket?" << std::endl;
            }
        }
        return false;
    }

    inline bool rollBackDisplacements(Key k, Data d, BFSQueue& bq)
    {
        hist[bq.size() -1] += 1;
        
        Key       k1;
        int       prev1;
        Bucket_t* b1;
        std::tie(k1, prev1, b1) = bq[bq.size()-1];

        Key       k2;
        int       prev2;
        Bucket_t* b2;
        while (prev1 >= 0)
        {
            std::tie(k2,prev2,b2) = bq[prev1];
            
            auto pop = b2->pop(k1);
            if (!pop.first)
            {
                std::cout << "serious issue with rollBack " << k1
                          << " pop " << pop.first << " " << pop.second << std::endl;
                for (size_t i = 0; i < tab.bs; ++i)
                {
                    std::cout << b2->p[i].first << " " << b2->p[i].second << std::endl;
                }
                return false; }
            if (!b1->insert(k1, pop.second))
            {   std::cout << "even more serious issue with rollBack" << std::endl; return false; }

            k1 = k2; prev1 = prev2; b1 = b2;
        }
        
        if (!b1->insert(k, d))
        {   std::cout << "failed final insert" << std::endl; return false; }
        
        return true;
        
            //return true;
    }
    
    inline bool insert(Key k, Data d, HashSplitter hash)
    {
        BFSQueue bq;
        Bucket_t* b1 = tab.getBucket1(hash);
        Bucket_t* b2 = tab.getBucket2(hash);

        bq.push_back(std::tuple<Key, int, Bucket_t*>(k, -1, b1));
        bq.push_back(std::tuple<Key, int, Bucket_t*>(k, -1, b2));
        
        for (size_t i = 0; i < visit; ++i)
        {
            if (expand(bq, i)) return rollBackDisplacements(k, d, bq);
        }
        
        return false;
    }
};

template<class Parent>
class dstrat_rwalk
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using Parent_t     = typename Parent::This_t;
    using HashSplitter = typename Parent::HashSplitter;
    using Bucket_t     = typename Parent::Bucket_t;

    Parent&      tab;
    std::mt19937 re;
    const size_t max_steps;
    size_t*      hist;

    dstrat_rwalk(Parent& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re((seed) ? seed : 30982391937209388ull), max_steps((steps) ? steps : 256)
    {
        hist = new size_t[max_steps];
        for (size_t i = 0; i < max_steps; ++i)
        {   hist[i] = 0;   }
    }
    ~dstrat_rwalk() { delete[] hist; }
    
    bool insert(Key k, Data d, HashSplitter hash)
    {        
        std::vector<std::tuple<Key, Data, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);

        auto    tk = k;
        auto    td = d;
        Bucket_t* tb;
        
        if (bin(re)) tb = tab.getBucket1(hash);
        else         tb = tab.getBucket2(hash);

        queue.emplace_back(tk,td,tb);

        for (size_t i = 0; tb->p[tab.bs-1].first && i<max_steps; ++i)
        {
            auto r = bsd(re);
            tk = tb->p[r].first;
            td = tb->p[r].second;
            auto hash = tab.h(tk);
            if (tab.getBucket1(hash) == tb) tb = tab.getBucket2(hash);
            else                          tb = tab.getBucket1(hash);

            queue.push_back(std::tuple<Key, Data, Bucket_t*>(tk,td,tb));


            // Explicit Cycle Detection (they are popped from the displacement queue)
            for (size_t j = 0; j < queue.size() - 1; ++j)
            {
                if (std::get<2>(queue[j]) == tb)
                {
                    while (queue.size() > j+1) queue.pop_back();
                    break;
                }
            }            
        }
        
        if (tb->p[tab.bs-1].first) { return false; }

        hist[queue.size() -1] += 1;
        
        for (size_t i = queue.size()-1; i > 0; --i)
        {            
            std::tie(tk,td,tb) = queue[i];
            if (! std::get<2>(queue[i-1])->remove(tk)) { std::cout << "e2" << std::endl; return false; }
            if (! tb->insert(tk, td)) { std::cout << "e1" << std::endl; return false; }           
        }

        if (! std::get<2>(queue[0])->insert(k, d)) { std::cout << "e3" << std::endl; return false; }
        return true;
    }
};


template<class Parent>
class dstrat_rwalk_cyclic
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using Parent_t     = typename Parent::This_t;
    using HashSplitter = typename Parent::HashSplitter;
    using Bucket_t     = typename Parent::Bucket_t;

    Parent&      tab;
    std::mt19937 re;
    const size_t max_steps;
    size_t*      hist;

    dstrat_rwalk_cyclic(Parent& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re((seed) ? seed : 30982391937209388ull), max_steps((steps) ? steps : 256)
    {
        hist = new size_t[max_steps];
        for (size_t i = 0; i < max_steps; ++i)
        {   hist[i] = 0;   }
    }
    ~dstrat_rwalk_cyclic() { delete[] hist; }
    
    bool insert(Key k, Data d, HashSplitter hash)
    {        
        std::vector<std::pair<std::pair<Key, Data>, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);

        auto tp = std::pair<Key,Data>(k, d);
        Bucket_t* tb;
        
        if (bin(re)) tb = tab.getBucket1(hash);
        else         tb = tab.getBucket2(hash);

        queue.emplace_back(tp,tb);
        auto ttp = tp;
        for (size_t i = 0; tb->p[tab.bs-1].first && i<max_steps; ++i)
        {
            auto r = bsd(re);
            tp = tb->p[r];
            tb->p[r] = ttp;
            ttp = tp;

            auto hash = tab.h(tp.first);
            if (tab.getBucket1(hash) == tb) tb = tab.getBucket2(hash);
            else                            tb = tab.getBucket1(hash);

            queue.emplace_back(tp,tb);
        }
        
        if (tb->insert(ttp.first, ttp.second))
        {
            hist[queue.size() -1] += 1;
            return true;
        }


        // No spot found -> revert all changes
        std::tie(ttp, tb) = queue[queue.size() - 1];
        for (size_t i = queue.size() - 2; i >= 1; --i)
        {
            std::tie(tp, tb) = queue[i];
            if (!(tb->remove(tp.first)))
            { std::cout << "f1" << std::endl; }
            if (!(tb->insert(ttp.first, ttp.second)))
            { std::cout << "f2" << std::endl; };
            ttp = tp;
        }
        
        return false;
    }
};


template<class, class, size_t>
class Bucket;


template<class K, class D, class H = std::hash<K>,
         template<class> class DS = dstrat_rwalk_cyclic,
         size_t TL = 128, size_t BS = 4>   // ALPHA (SIZE CONSTRAINT COULD BE PARAMETER BUT INTEGRAL TYPE PROBLEM)
class SpaceGrow
{
public:
    using Key       = K;
    using Data      = D;
    using FRet      = std::pair<bool, Data>;    

    SpaceGrow(size_t _capacity = 0, double size_constraint = 1.1, size_t dis_steps = 0, size_t seed = 0)
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

    ~SpaceGrow() = default;

    SpaceGrow(const SpaceGrow&) = delete;
    SpaceGrow& operator=(const SpaceGrow&) = delete;

    SpaceGrow(SpaceGrow&& rhs) = default;
    SpaceGrow& operator=(SpaceGrow&& rhs) = default;

    bool insert(Key k, Data d);
    FRet find  (Key k);
    bool remove(Key k);

    void printHist();
    
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
    size_t       curGrowTable;;
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

    Bucket() { for (size_t i = 0; i < BS; ++i) p[i] = std::make_pair(Key(), Data());}
    Bucket(const Bucket& rhs) = default;
    Bucket& operator=(const Bucket& rhs) = default;
    
    bool   insert(Key k, Data d);
    FRet   find  (Key k);
    bool   remove(Key k);
    FRet   pop   (Key k);
    
    int    probe (Key k);
    
    std::pair<Key, Data> p[BS];
};




template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::insert(Key k, Data d)
{
    auto hash = h(k);
        
    auto p1 = getBucket1(hash)->probe(k);//llt[ hash.tab1 ].probe( k, hash.loc1 );
    auto p2 = getBucket2(hash)->probe(k);//llt[ hash.tab2 ].probe( k, hash.loc2 );
    
    if ((p1 < 0) || (p2 < 0)) return false;

    auto r = false;
    
    // SHOULD CHECK, IF ALREADY INCLUDED
    if (p1 > p2)
    {
        // insert into p1
        r = getBucket1(hash)->insert(k,d);
    }
    else if (p2 > p1)
    {
        // insert into p2
        r = getBucket2(hash)->insert(k,d);
    }
    else if (p1)
    {
        // insert into p1
        r = getBucket1(hash)->insert(k,d);
    }
    else
    {
        // no space => displace stuff
        r = displacer.insert(k, d, hash);
    }
    
    if (r) incElements();

    return r;
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
typename SpaceGrow<K,D,HF,DS,TL,BS>::FRet SpaceGrow<K,D,HF,DS,TL,BS>::find(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->find(k); //llt[hash.tab1].find(k, hash.loc1);

    if (p1.first) return p1;
    else          return getBucket2(hash)->find(k); //llt[hash.tab2].find(k, hash.loc2);
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
bool SpaceGrow<K,D,HF,DS,TL,BS>::remove(Key k)
{
    auto hash = h(k);
    auto p1 = getBucket1(hash)->remove(k); //llt[hash.tab1].remove(k, hash.loc1);
    
    if (p1) return true;
    else    return getBucket2(hash)->remove(k); //llt[hash.tab2].remove(k,hash.loc2);
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
            auto e    = curr->p[j];
            if (! e.first) break;
            auto hash = h(e.first);
            if      (getBucket1(hash) == curr) target[bitmask & hash.loc1].insert(e.first, e.second);
            else if (getBucket2(hash) == curr) target[bitmask & hash.loc2].insert(e.first, e.second);
            else { std::cout << "something is wrong neither in first, nor second bucket." << std::endl; }
        }
    }
}

template<class K, class D, class HF, template<class> class DS, size_t TL, size_t BS>
void SpaceGrow<K,D,HF,DS,TL,BS>::printHist()
{
    size_t gHist[BS];
    for (size_t i = 0; i < BS; ++i) gHist[i] = 0;

    out ("# tab", 5);
    out ("full" , 6);
    out ("1 sp" , 6);
    out ("2 sp" , 6);
    out ("empt" , 6);
    out ("n"    , 8);
    out ("cap"  , 8);
    std::cout << std::endl;
    
    for (size_t tl = 0; tl < TL; ++tl)
    {
        size_t lHist[BS];
        for (size_t i = 0; i < BS; ++i) lHist[i] = 0;
        
        for (size_t j = 0; j <= llb[tl]; ++j)
        {
            auto a = llt[tl][j].probe(1);
            if (a >= 0) ++lHist[a];
        }

        size_t n = 0;
        out (tl, 5);
        for (size_t i = 0; i < BS; ++i)
        {
            out (lHist[i], 6);
            n += lHist[i] * (BS - i);
            gHist[i] += lHist[i];
        }
        out (n, 8);
        out (llb[tl]*BS, 8);
        std::cout << std::endl;
    }

    size_t n = 0;
    out ("#all", 5);
    for (size_t i = 0; i < BS; ++i)
    {
        out (gHist[i], 6);
        n += gHist[i] * (BS - i);
    }
    out (n, 8);
    out (curGrowAmount * (TL+curGrowTable), 8);
    std::cout << std::endl;    
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
        if ( p[i].first == k ) return std::make_pair(true, p[i].second);
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
        if (p[i].first == k)
        {
            Data d = p[i].second;
            for (size_t j = i+1; j < BS; ++j)
            {
                if (p[j].first) { p[i] = p[j]; i = j; }
                else break;
            }
            p[i] = std::make_pair(Key(), Data());
            return std::make_pair(true, d);                                    
        }
        else if (! p[i].first)
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
        if (!p[i].first)      return BS - i;
        if ( p[i].first == k) return -1;
        
    }
    return 0;
}
