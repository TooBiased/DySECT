//#include "include/spacegrow.h"
#include "selection.h"

#include "include/strategies/dstrat_bfs.h"
#include "include/strategies/dstrat_rwalk.h"
#include "include/strategies/dstrat_rwalk_cyclic.h"

#include "utils/hashfct.h"
#include "utils/commandline.h"

#include <random>
#include <iostream>
#include <fstream>

template <class T>
inline void print(std::ostream& out, const T& t, size_t w)
{
    out.width(w);
    out << t << " " << std::flush;
}

template<class Table>
void print_hist(Table& table, std::string name)
{
    std::ofstream out(name, std::ofstream::out);

    auto& hcount = table.hcounter;

    print(out, "# steps", 7);
    print(out, "nFitted", 8);
    out << std::endl;

    for (size_t i = 0; i < hcount.steps; ++i)
    {
        print(out, i, 7);
        print(out, hcount.hist[i], 8);
        out << std::endl;
    }
    out.close();
}


template<class Table>
void print_dist(Table& table, std::string name)
{
    std::ofstream out(name, std::ofstream::out);

    print (out, "# tab", 5);
    size_t gHist[table.bs+1];
    for (size_t i = 0; i <= table.bs; ++i)
    {
        gHist[i] = 0;
        print (out, i, 6);
    }

    print (out, "n"    , 8);
    print (out, "cap"  , 8);
    out << std::endl;

    for (size_t tl = 0; tl < table.tl; ++tl)
    {
        size_t lHist[table.bs+1];
        for (size_t i = 0; i <= table.bs; ++i) lHist[i] = 0;

        auto ltab = table.getTable(tl);

        for (size_t j = 0; j < ltab.first; ++j)
        {
            auto a = ltab.second[j].probe(0);
            if (a >= 0) ++lHist[a];
        }

        size_t n = 0;
        print (out, tl, 5);
        for (size_t i = 0; i <= table.bs; ++i)
        {
            print (out, lHist[i], 6);
            n += lHist[i] * (table.bs - i);
            gHist[i] += lHist[i];
        }
        print (out, n, 8);
        print (out, ltab.first*table.bs, 8);
        out << std::endl;
    }

    size_t n = 0;
    print (out, "#all", 5);
    for (size_t i = 0; i <= table.bs; ++i)
    {
        print (out, gHist[i], 6);
        n += gHist[i] * (table.bs - i);
    }
    print (out, n, 8);
    //print (out, table.curGrowAmount * (table.tl+table.curGrowTable), 8);
    out << std::endl;

    out.close();
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

    auto errors = 0;
    bool first  = true;

    HASHTYPE<size_t, size_t, HASHFCT, Displacer, TL, BS> table(cap, alpha, steps);

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

    print_dist(table, name + ".dist");
    print_hist(table, name + ".hist");

    delete[] keys;

    return 0;
}

template<template<class> class Displacer, size_t TL>
int test (size_t n, size_t cap, size_t steps, double alpha, std::string name, size_t bs)
{
    switch (bs)
    {
        //  case 4:
        //return test<Displacer, TL, 4 >(n,cap,steps,alpha,name);
    //case 6:
        //return test<Displacer, TL, 6>(n,cap,steps,alpha,name);
    case 8:
        return test<Displacer, TL, 8 >(n,cap,steps,alpha,name);
        //case 16:
        //return test<Displacer, TL, 16>(n,cap,steps,alpha,name);
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
    case 2:
        return test<Displacer, 2>(n,cap,steps,alpha,name,bs);

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
        //case 64:
        //return test<Displacer, 64 >(n,cap,steps,alpha,name,bs);
    case 128:
        return test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        //case 256:
        //return test<Displacer, 256>(n,cap,steps,alpha,name,bs);
        /*
    case 512:
        test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        break;*/
        //case 2048:
        //  return test<Displacer, 2048>(n,cap,steps,alpha,name,bs);

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
