#include "include/spacegrow.h"
#include "include/strategies/dstrat_bfs.h"
#include "include/strategies/dstrat_rwalk.h"
#include "include/strategies/dstrat_rwalk_cyclic.h"

#include "utils/hashfct.h"
#include "utils/commandline.h"

#include <random>
#include <iostream>
#include <fstream>



template<template<class> class Displacer>
int test(size_t n, size_t cap, size_t steps, double alpha, std::string name)
{    
    constexpr size_t range = (1ull<<63) -1;

    size_t* keys = new size_t[n];

    std::uniform_int_distribution<uint64_t> dis(1,range);
    std::mt19937_64 re;
    
    for (size_t i = 0; i < n; ++i)
    {
        keys[i] = dis(re);
    }
    
    std::cout << "randomness generated" << std::endl;

    auto errors = 0;
    bool first  = true;
    
    SpaceGrow<size_t, size_t, murmur_hasher, Displacer> table(cap, alpha, steps);
    
    std::cout << "table generated"      << std::endl;
    
    for (size_t i = 0; i < n; ++i)
    {
        if (!table.insert(keys[i], i))
        {
            ++errors;
            
            if (first)
            {
                first = false;
                std::cout << "first error at i=" << i << std::endl;
            }
        }
    }

    std::cout << "inserted elements encountered " << errors << " errors" << std::endl;
    
    auto count = 0;
    errors = 0;
    for (size_t i = 0; i < n; ++i)
    {
        auto e = table.find(keys[i]);
        if (e.first && (e.second == i)) count++;
        else if (e.first) std::cout << "wat" << std::endl;
        else errors++;
    }

    std::cout << "count: " << count << "  errors: " << errors << std::endl;

    
    std::ofstream distOut(name + ".dist", std::ofstream::out);
    table.printDist(distOut);
    distOut.close();

    std::ofstream histOut(name + ".hist", std::ofstream::out);
    table.printHist(histOut);
    histOut.close();    

    delete[] keys;
    
    return 0;
}



int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      n     = c.intArg("-n"    , 1000000);
    const size_t      cap   = c.intArg("-cap"  , 0);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "temp");
    const double      alpha = c.doubleArg("-alpha", 1.1);

    if      (c.boolArg("-bfs"))
    {
        return test<dstrat_bfs>         (n, cap, steps, alpha, name);        
    }
    else if (c.boolArg("-rwalk"))
    {
        return test<dstrat_rwalk>       (n, cap, steps, alpha, name);        
    }
    else if (c.boolArg("-rwalkcyc"))
    {
        return test<dstrat_rwalk_cyclic>(n, cap, steps, alpha, name);        
    }

    std::cout << "no displacement strategy chosen" << std::endl;
    test<dstrat_triv>(n, cap, steps, alpha, name);

    return 1;
}
