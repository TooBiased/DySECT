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
    using Table = HASHTYPE<size_t, size_t, HASHFCT, Config>;

    inline void print_headline(std::ostream& out)
    {
        print(out, "it"    , 3);
        print(out, "alpha" , 5);
        print(out, "bsize" , 5);
        print(out, "ntabl" , 5);
        print(out, "cap"   , 8);
        print(out, "n_pre" , 8);
        print(out, "n_main", 8);
        print(out, "t_pre" , 8);
        print(out, "t_in"  , 8);
        print(out, "t_find", 8);
        print(out, "unsucc", 6);
        print(out, "lost"  , 6);
        out << std::endl;
    }

    inline void print_timing(std::ostream& out, size_t i,
                             double alpha, size_t bs  , size_t tl,
                             size_t cap  , size_t pre , size_t n ,
                             double d_pre, double d_in, double d_fn,
                             size_t unsucc, size_t lost_elem)
    {
        print(out, i        , 3);
        print(out, alpha    , 5);
        print(out, bs       , 5);
        print(out, tl       , 5);
        print(out, cap      , 8);
        print(out, pre      , 8);
        print(out, n        , 8);
        print(out, d_pre    , 8);
        print(out, d_in     , 8);
        print(out, d_fn     , 8);
        print(out, unsucc   , 6);
        print(out, lost_elem, 6);
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t pre,  size_t cap, size_t steps, double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n+pre];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n+pre; ++i)
        {
            keys[i] = dis(re);
        }

        std::ofstream file(name + ".time", std::ofstream::out | std::ofstream::app);
        print_headline(file);

        for (size_t i = 0; i < it; ++i)
        {
            auto errors = 0;

            Table table(cap, alpha, steps);

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < pre; ++i)
            {
                table.insert(keys[i], i);
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            for (size_t i = pre; i < pre + n; ++i)
            {
                if (!table.insert(keys[i], i))
                {
                    ++errors;
                }
            }

            auto t2 = std::chrono::high_resolution_clock::now();
            auto count = 0;
            for (size_t i = pre; i < pre + n; ++i)
            {
                auto e = table.find(keys[i]);
                if (e.first && (e.second == i)) count++;
            }
            auto t3 = std::chrono::high_resolution_clock::now();

            double d_pre = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1000.;
            double d_in  = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()/1000.;
            double d_fn  = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count()/1000.;

            print_timing(file, it, alpha, Config::bs, Config::tl, cap, pre, n,
                         d_pre, d_in, d_fn, errors, n - errors -count);
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      it    = c.intArg("-it"   , 5);
    const size_t      pre   = c.intArg("-pre"  , 1000000);
    const size_t      n     = c.intArg("-n"    , 1000000);
    const size_t      cap   = c.intArg("-cap"  , pre + n);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "temp");
    const double      alpha = c.doubleArg("-alpha", 1.1);

    return Chooser::execute<Test,no_hist_count> (c, it, n, pre, cap, steps, alpha, name);
}
