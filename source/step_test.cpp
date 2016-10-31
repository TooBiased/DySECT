#include "include/spacegrow.h"
#include "include/strategies/dstrat_bfs.h"
#include "include/strategies/dstrat_rwalk.h"
#include "include/strategies/dstrat_rwalk_cyclic.h"

#include "utils/hashfct.h"
#include "utils/commandline.h"

#include <random>
#include <iostream>
#include <fstream>

template<class Table>
size_t catchup(Table& table, size_t* keys, size_t i)
{
    for (size_t i = 0; i < n; ++i)
    {
        table.insert(keys[i], i);
    }
}


template<template<class> class Displacer, size_t TL, size_t BS>
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

    SpaceGrow<size_t, size_t, HASHFCT, Displacer, TL, BS> table_d(cap, alpha, steps);
    SpaceGrow<size_t, size_t, HASHFCT, Displacer, TL, BS> table_s(table_d.capacity, 1.0, steps);

    auto curr_d_size = table_d.capacity;
    auto errors_d = 0;
    auto errors_s = 0;

    std::cout << "tables generated"      << std::endl;

    for (size_t i = 0; i < n; ++i)
    {
        if (!table_d.insert(keys[i], i)) { ++errors_d; }
        if (!table_s.insert(keys[i], i)) { ++errors_s; }
        if ( table_d.capacity != curr_d_size )
        {
            print_stuff(table_d, table_s, );
            table_s = std::move(table_s(table_d.capacity, 1.0, steps));
            catchup(table_s, keys, i);
            table_s.clearHist();

        }
    }

    std::cout << "inserted elements encountered " << errors << " errors" << std::endl;


    std::ofstream distOut(name + ".dist", std::ofstream::out);
    table.printDist(distOut);
    distOut.close();

    std::ofstream histOut(name + ".hist", std::ofstream::out);
    table.printHist(histOut);
    histOut << steps+100 << "     " << errors << std::endl;
    histOut.close();


    delete[] keys;

    return 0;
}

template<template<class> class Displacer, size_t TL>
int test (size_t n, size_t cap, size_t steps, double alpha, std::string name, size_t bs)
{
    switch (bs)
    {
    case 4:
        return test<Displacer, TL, 4 >(n,cap,steps,alpha,name);
    //case 6:
        //return test<Displacer, TL, 6>(n,cap,steps,alpha,name);
    case 8:
        return test<Displacer, TL, 8 >(n,cap,steps,alpha,name);
    case 16:
        return test<Displacer, TL, 16>(n,cap,steps,alpha,name);
    default:
        std::cout << "UNKNOWN BS " << bs << std::endl;
        return 32;
    }
}

template<template<class> class Displacer>
int test (size_t n, size_t cap, size_t steps, double alpha, std::string name, size_t tl, size_t bs)
{
    switch (tl)
    {
        /*
    case 8:
        test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        break;
    case 16:
        test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        break;
    case 32:
        test<Displacer, 32 >(n,cap,steps,alpha,name,bs);
        break;
        */
    case 64:
        return test<Displacer, 64 >(n,cap,steps,alpha,name,bs);
    case 128:
        return test<Displacer, 128>(n,cap,steps,alpha,name,bs);
    case 256:
        return test<Displacer, 256>(n,cap,steps,alpha,name,bs);
        /*
    case 512:
        test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        break;*/
    case 2048:
        return test<Displacer, 2048>(n,cap,steps,alpha,name,bs);

    default:
        std::cout << "UNKNOWN TL " << tl << std::endl;
        return 16;
    }
}

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      n     = c.intArg("-n"    , 1000000);
    const size_t      cap   = c.intArg("-cap"  , 0);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "temp");
    const double      alpha = c.doubleArg("-alpha", 1.1);
    const size_t      tl    = c.intArg("-tl"   , 128);
    const size_t      bs    = c.intArg("-bs"   , 4);

    if      (c.boolArg("-bfs"))
    {
        return test<dstrat_bfs>         (n, cap, steps, alpha, name, tl, bs);
    }
    else if (c.boolArg("-rwalk"))
    {
        return test<dstrat_rwalk>       (n, cap, steps, alpha, name, tl, bs);
    }
    else if (c.boolArg("-rwalkcyc"))
    {
        return test<dstrat_rwalk_cyclic>(n, cap, steps, alpha, name, tl, bs);
    }

    std::cout << "ERROR: choose Displacement Strategy" << std::endl;

    return 1;
}
