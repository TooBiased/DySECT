#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>



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
    const size_t steps;
    size_t*      hist;

    dstrat_rwalk(Parent& parent, size_t steps=256, size_t seed=30982391937209388ull)
        : tab(parent), re((seed) ? seed : 30982391937209388ull), steps((steps) ? steps : 256)
    {
        hist = new size_t[steps];
        for (size_t i = 0; i < steps; ++i)
        {   hist[i] = 0;   }
    }
    
    ~dstrat_rwalk() { delete[] hist; }
    
    bool insert(std::pair<Key,Data> t, HashSplitter hash)
    {        
        std::vector<std::pair<std::pair<Key, Data>, Bucket_t*> > queue;
        std::uniform_int_distribution<size_t> bin(0,1);
        std::uniform_int_distribution<size_t> bsd(0,tab.bs-1);

        auto      tp = t;
        Bucket_t* tb;
        
        if (bin(re)) tb = tab.getBucket1(hash);
        else         tb = tab.getBucket2(hash);

        queue.emplace_back(tp,tb);

        for (size_t i = 0; !(tb->space()) && i<steps; ++i)
        {
            auto r = bsd(re);
            tp = tb->get(r);
            auto hash = tab.h(tp.first);
            if (tab.getBucket1(hash) == tb) tb = tab.getBucket2(hash);
            else                            tb = tab.getBucket1(hash);

            queue.emplace_back(tp,tb);


            // Explicit Cycle Detection (they are popped from the displacement queue)
            for (size_t j = 0; j < queue.size() - 1; ++j)
            {
                if (queue[j].second == tb)
                {
                    while (queue.size() > j+1) queue.pop_back();
                    break;
                }
            }            
        }
        
        if (! tb->space()) { return false; }

        hist[queue.size() -1] += 1;
        
        for (size_t i = queue.size()-1; i > 0; --i)
        {            
            std::tie(tp,tb) = queue[i];
            if (! queue[i-1].second->remove(tp.first)) { std::cout << "e2" << std::endl; return false; }
            if (! tb->insert(tp.first, tp.second))      { std::cout << "e1" << std::endl; return false; }           
        }

        if (! queue[0].second->insert(t)) { std::cout << "e3" << std::endl; return false; }
        return true;
    }
};
