#include "include/spacegrow.h"
#include "utils/hashfct.h"
#include "utils/commandline.h"

#include <random>
#include <iostream>


constexpr size_t Nbr   = 10000000;
constexpr size_t Alpha = 1.1;
constexpr size_t Steps = 1024;
constexpr size_t range = (1ull<<63) -1;



template<typename T>
inline void out(const T& t, int space)
{
    std::cout.width(space);
    std::cout << t << " " << std::flush;
}



int main(int argn, char** argc)
{
    CommandLine c{argn, argc};
    size_t n     = c.intArg   ("-n"    , Nbr  );
    double alpha = c.doubleArg("-a"    , Alpha);
    size_t cap   = c.intArg   ("-cap"  , Nbr);
    size_t steps = c.intArg   ("-steps", Steps);

    size_t* keys = new size_t[n];

    std::uniform_int_distribution<uint64_t> dis(1,range);
    std::mt19937_64 re;
    
    for (size_t i = 0; i < n; ++i)
    {
        keys[i] = dis(re);
    }

    size_t errors[12]; for (size_t i = 0; i<12; ++i) errors[i] = 0;
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 256, 4> tab0 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 64 , 4> tab1 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 16 , 4> tab2 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 256, 8> tab3 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 64 , 8> tab4 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_bfs  , 16 , 8> tab5 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 256, 4> tab6 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 64 , 4> tab7 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 16 , 4> tab8 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 256, 8> tab9 (cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 64 , 8> tab10(cap, alpha, steps);
    SpaceGrow<size_t, size_t, murmur_hasher, dstrat_rwalk, 16 , 8> tab11(cap, alpha, steps);

    for (size_t i = 0; i < n; ++i)
    {
        auto key = keys[i];
        
        if ( !tab0 .insert(key, i) ) ++errors [0];
        if ( !tab1 .insert(key, i) ) ++errors [1];
        if ( !tab2 .insert(key, i) ) ++errors [2];
        if ( !tab3 .insert(key, i) ) ++errors [3];
        if ( !tab4 .insert(key, i) ) ++errors [4];
        if ( !tab5 .insert(key, i) ) ++errors [5];
        if ( !tab6 .insert(key, i) ) ++errors [6];
        if ( !tab7 .insert(key, i) ) ++errors [7];
        if ( !tab8 .insert(key, i) ) ++errors [8];
        if ( !tab9 .insert(key, i) ) ++errors [9];
        if ( !tab10.insert(key, i) ) ++errors[10];
        if ( !tab11.insert(key, i) ) ++errors[11];
    }

    // table.printHist();

    for (size_t i = 0; i < steps; ++i)
    {
        out(i, 5);
        out(tab0 .displacer.hist[i], 7);
        out(tab1 .displacer.hist[i], 7);
        out(tab2 .displacer.hist[i], 7);
        out(tab3 .displacer.hist[i], 7);
        out(tab4 .displacer.hist[i], 7);
        out(tab5 .displacer.hist[i], 7);
        out(tab6 .displacer.hist[i], 7);
        out(tab7 .displacer.hist[i], 7);
        out(tab8 .displacer.hist[i], 7);
        out(tab9 .displacer.hist[i], 7);
        out(tab10.displacer.hist[i], 7);
        out(tab11.displacer.hist[i], 7);
        std::endl;
    }

    delete[] keys;
    
    

    return 0;
}
