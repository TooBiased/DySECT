//#include "include/spacegrow.h"
#include "selection.h"

// #include "include/strategies/dstrat_bfs.h"
// #include "include/strategies/dstrat_rwalk.h"
// #include "include/strategies/dstrat_rwalk_cyclic.h"

#include "utils/hashfct.h"
#include "utils/commandline.h"

#ifdef MALLOC_COUNT
#include "malloc_count.h"
#endif

#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

template <class T>
inline void print(std::ostream& out, const T& t, size_t w)
{
    out.width(w);
    out << t << " " << std::flush;
}

template<class Config>
struct Test
{
    //using Table = ProbIndependentBase<HASHTYPE<size_t, size_t, HASHFCT, Config> >;
    using Table = HASHTYPE<size_t, size_t, HASHFCT, Config>;

    inline void print_headline(std::ostream& out)
    {
        print(out, "# it"   , 4);
        print(out, "alpha"  , 5);
        Table::print_init_header(out);
        print(out, "cap"    , 9);
        print(out, "n_0"  , 9);
        print(out, "n_full" , 9);
        print(out, "t_in0"  , 8);
        print(out, "t_in1"  , 8);
        print(out, "t_find+", 8);
        print(out, "t_find-", 8);
        print(out, "in_err" , 6);
        print(out, "fi_err" , 6);
        #ifdef MALLOC_COUNT
        print(out, "memory" , 7);
        #endif
        out << std::endl;
    }

    inline void print_timing(std::ostream& out,
                             size_t i,
                             double alpha,
                             size_t cap  ,
                             size_t n0 ,
                             size_t n ,
                             double d_in0,
                             double d_in1,
                             double d_fn0,
                             double d_fn1,
                             size_t in_err,
                             size_t fi_err,
                             Table& table)
    {
        print(out, i     , 4);
        print(out, alpha , 5);
        table.print_init_data(out);
        print(out, cap   , 9);
        print(out, n0    , 9);
        print(out, n     , 9);
        print(out, d_in0 , 8);
        print(out, d_in1 , 8);
        print(out, d_fn0 , 8);
        print(out, d_fn1 , 8);
        print(out, in_err, 6);
        print(out, fi_err, 6);
        #ifdef MALLOC_COUNT
        print(out, double(malloc_count_current())/double(8*2*n)-1., 7);
        #endif
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t n0,  size_t cap, size_t steps,
                   double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[2*n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < 2*n; ++i)
        {
            keys[i] = dis(re);
        }

        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".time",
                                      std::ofstream::out | std::ofstream::app);
        //std::ofstream file(name + ".time",
        //                   std::ofstream::out | std::ofstream::app);
        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            Table table(cap, alpha, steps);

            auto in_errors = 0ull;
            auto fin_errors = 0ull;


            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < n0; ++i)
            {
                if (!table.insert(keys[i], i).second) ++in_errors;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            for (size_t i = n0; i < n; ++i)
            {
                if (!table.insert(keys[i], i).second) ++in_errors;
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            //const Table& ctable = table;
            for (size_t i = 0; i < n; ++i)
            {
                auto e = table.find(keys[i]);
                if ( (e == table.end()) || ((*e).second != i))
                    //if (ctable.at(keys[i]) != i)
                    fin_errors++;
            }
            auto t3 = std::chrono::high_resolution_clock::now();
            for (size_t i = n; i < 2*n; ++i)
            {
                auto e = table.find(keys[i]);
                if ( (e != table.end()) && (keys[(*e).second] != keys[i]))
                    fin_errors++;
            }
            auto t4 = std::chrono::high_resolution_clock::now();

            double d_in0 = std::chrono::duration_cast<std::chrono::microseconds> (t1 - t0).count()/1000.;
            double d_in1 = std::chrono::duration_cast<std::chrono::microseconds> (t2 - t1).count()/1000.;
            double d_fn0 = std::chrono::duration_cast<std::chrono::microseconds> (t3 - t2).count()/1000.;
            double d_fn1 = std::chrono::duration_cast<std::chrono::microseconds> (t4 - t3).count()/1000.;

            print_timing(*file, i, alpha, cap, n0, n, d_in0, d_in1, d_fn0, d_fn1,
                         in_errors, fin_errors, table);

            // Iterator Test
            // size_t temp0 = 0;
            // size_t temp1 = 0;
            // for (auto it = table.begin(); it != table.end(); it++)
            // {
            //     if (it->second < 100000)
            //         ++temp1;
            //     ++temp0;
            // }
            // std::cout << temp1 << "/" << temp0 << std::endl;
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      it    = c.intArg("-it"   , 5);
    const size_t      n     = c.intArg("-n"    , 2000000);
    const size_t      n0    = c.intArg("-pre"  , n/2);
    if (n0 > n) std::cout << "n0 (pre) has to be smaller than n! set n = n0 = " << n0 << std::endl;
    const size_t      cap   = c.intArg("-cap"  , n);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "");
    const double      alpha = c.doubleArg("-alpha", 1.1);
    const double      load  = c.doubleArg("-load" , 2.0);
    const double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,no_hist_count> (c, it, n, n0, cap, steps, alpha, name);
}
