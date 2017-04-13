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
        print(out, "pat"    , 4);
        print(out, "pre"    , 9);
        print(out, "n"      , 9);
        print(out, "t_pre"  , 8);
        print(out, "t_mix"  , 8);
        print(out, "in_err" , 6);
        #ifdef MALLOC_COUNT
        print(out, "memory" , 7);
        #endif
        out << std::endl;
    }

    inline void print_timing(std::ostream& out,
                             size_t i,
                             double alpha,
                             size_t cap  ,
                             size_t pattern,
                             size_t pre ,
                             size_t n ,
                             double d_pre,
                             double d_mix,
                             size_t errors,
                             Table& table)
    {
        print(out, i     , 4);
        print(out, alpha , 5);
        table.print_init_data(out);
        print(out, cap   , 9);
        print(out, pattern, 4);
        print(out, pre    , 9);
        print(out, n     , 9);
        print(out, d_pre , 8);
        print(out, d_mix , 8);
        print(out, errors, 6);
        #ifdef MALLOC_COUNT
        print(out, double(malloc_count_current()-(n+pre)*8)/
                   (0.1*double(pattern*8*2*n)+double(8*2*pre)), 7);
        #endif
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t pre, size_t cap, size_t pattern, size_t steps,
                   double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[pre+n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < pre; ++i)
        {
            keys[i] = dis(re);
        }
        size_t pcount = 0;
        for (size_t i = pre; i < n+pre; ++i)
        {
            if (pcount < pattern)
            {
                keys[i] = dis(re);
            }
            else
            {
                std::uniform_int_distribution<uint64_t> find_dis(0,i-1);
                keys[i] = keys[find_dis(re)];
            }
            pcount  = (pcount < 9) ? pcount+1 : 0;
        }


        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".mix",
                                      std::ofstream::out | std::ofstream::app);

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            Table table(cap, alpha, steps);

            auto errors  = 0ull;
            auto ferrors = 0ull;

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < pre; ++i)
            {
                if (!table.insert(keys[i], i).second) ++errors;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            pcount = 0;
            for (size_t i = pre; i < n+pre; ++i)
            {
                if (pcount < pattern)
                {
                    if (!table.insert(keys[i], i).second) ++errors;
                }
                else
                {
                    auto res = table.find(keys[i]);
                    if (res == table.end()) ++ferrors;
                }
                pcount  = (pcount < 9) ? pcount+1 : 0;
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            double d_pre = std::chrono::duration_cast<std::chrono::microseconds> (t1 - t0).count()/1000.;
            double d_mix = std::chrono::duration_cast<std::chrono::microseconds> (t2 - t1).count()/1000.;

            print_timing(*file, i, alpha, cap, pattern, pre, n, d_pre, d_mix, errors, table);
            if (ferrors && (!errors)) (*file) << "Critical Error" << std::endl;
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    size_t      it    = c.intArg("-it"   , 5);
    size_t      n     = c.intArg("-n"    , 1000000);
    size_t      n0    = c.intArg("-pre"  , n/4);
    size_t      cap   = c.intArg("-cap"  , n0);
    size_t      pattern = c.intArg("-pattern", 5);
    size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "");
    double      alpha = c.doubleArg("-alpha", 1.1);
    double      load  = c.doubleArg("-load" , 2.0);
    double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,no_hist_count> (c, it, n, n0, cap, pattern, steps, alpha, name);
}
