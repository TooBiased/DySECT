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
void print_hist(Table& table, size_t lvl, std::string name)
{
    std::ofstream out(name, std::ofstream::app);

    auto& disp = table.displacer;

    print(out, "# lvl", 5);
    print(out, "steps", 7);
    print(out, "nFitted", 8);
    out << std::endl;

    for (size_t i = 0; i < disp.steps; ++i)
    {
        print(out, lvl, 5);
        print(out, i, 7);
        print(out, disp.hist[i], 8);
        out << std::endl;
    }
    out.close();
}

template<class Table>
void print_dist(Table& table, size_t lvl, std::string name)
{
    std::ofstream out(name, std::ofstream::app);

    print (out, "# lvl", 5);
    print (out, "tab", 5);
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

        auto ctab = table.getTable(tl);

        for (size_t j = 0; j < ctab.first; ++j)
        {
            auto a = ctab.second[j].probe(0);
            if (a >= 0) ++lHist[a];
        }

        size_t n = 0;
        print (out, lvl, 5);
        print (out, tl, 5);
        for (size_t i = 0; i <= table.bs; ++i)
        {
            print (out, lHist[i], 6);
            n += lHist[i] * (table.bs - i);
            gHist[i] += lHist[i];
        }
        print (out, n, 8);
        print (out, ctab.first*table.bs, 8);
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

template<class Table>
void print_stuff(Table& tab_d, Table& tab_s, size_t lvl, std::string name)
{
    print_hist(tab_d, lvl, name+"_dynamic.hist");
    print_hist(tab_s, lvl, name+"_static.hist");
    print_dist(tab_d, lvl, name+"_dynamic.dist");
    print_dist(tab_s, lvl, name+"_static.dist");
}

void init_output_files(std::string name)
{
    std::ofstream out0(name+"_dynamic.hist", std::ofstream::out);
    out0.close();
    std::ofstream out1(name+"_static.hist" , std::ofstream::out);
    out1.close();
    std::ofstream out2(name+"_dynamic.dist", std::ofstream::out);
    out2.close();
    std::ofstream out3(name+"_static.dist" , std::ofstream::out);
    out3.close();
}

template<class Table>
void catchup(Table& table, size_t* keys, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        table.insert(keys[i], i);
    }
}


template<template<class> class Displacer, size_t TL, size_t BS>
int test(size_t n, size_t cap, size_t steps, double alpha, std::string name)
{
    init_output_files(name);

    constexpr size_t range = (1ull<<63) -1;

    size_t* keys = new size_t[n];

    std::uniform_int_distribution<uint64_t> dis(1,range);
    std::mt19937_64 re;

    for (size_t i = 0; i < n; ++i)
    {
        keys[i] = dis(re);
    }

    std::cout << "randomness generated" << std::endl;

    HASHTYPE<size_t, size_t, HASHFCT, Displacer, TL, BS> table_d(cap, alpha, steps);
    HASHTYPE<size_t, size_t, HASHFCT, Displacer, TL, BS> table_s(table_d.capacity, 1.0, steps);

    auto curr_d_size = table_d.capacity;
    auto errors_d = 0;
    auto errors_s = 0;
    auto lvl = 0;

    std::cout << "tables generated"      << std::endl;

    for (size_t i = 0; i < n; ++i)
    {
        if (!table_d.insert(keys[i], i)) { ++errors_d; }
        if (!table_s.insert(keys[i], i)) { ++errors_s; }
        if ( table_d.capacity != curr_d_size )
        {
            print_stuff(table_d, table_s, lvl++, name);
            curr_d_size = table_d.capacity;

            std::cout << "NEW TABLE " << lvl << " size:" << curr_d_size << std::endl;

            table_s = std::move(HASHTYPE<size_t, size_t, HASHFCT, Displacer, TL, BS>
                                (curr_d_size, 1.0, steps));
            catchup(table_s, keys, i);
            table_s.clearHist();
            table_d.clearHist();
        }
    }

    std::cout << "inserted elements encountered " << errors_d
              << " " << errors_s << " errors" << std::endl;

    delete[] keys;

    return 0;
}

template<template<class> class Displacer, size_t TL>
int test (size_t n, size_t cap, size_t steps, double alpha, std::string name, size_t bs)
{
    switch (bs)
    {
        //case 4:
        //return test<Displacer, TL, 4 >(n,cap,steps,alpha,name);
    //case 6:
        //return test<Displacer, TL, 6>(n,cap,steps,alpha,name);
    case 8:
        return test<Displacer, TL, 8 >(n,cap,steps,alpha,name);
        //case 16:
        // return test<Displacer, TL, 16>(n,cap,steps,alpha,name);
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
        //case 64:
        //   return test<Displacer, 64 >(n,cap,steps,alpha,name,bs);
    case 128:
        return test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        //case 256:
        // return test<Displacer, 256>(n,cap,steps,alpha,name,bs);
        /*
    case 512:
        test<Displacer, 128>(n,cap,steps,alpha,name,bs);
        break;*/
        //case 2048:
        //return test<Displacer, 2048>(n,cap,steps,alpha,name,bs);

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
    const size_t      bs    = c.intArg("-bs"   , 8);

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
