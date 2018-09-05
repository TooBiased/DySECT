//#include "include/spacegrow.h"
#include "selection.h"

#include "utils/default_hash.h"
#include "utils/commandline.h"
#include "utils/thread_basics.h"

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
    //using Table = ProbIndependentBase<HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config> >;
    using Table = HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config>;

    constexpr static size_t block_size = 100000;

    inline void print_headline(std::ostream& out)
    {
        print(out, "# it"   , 4);
        print(out, "block"  , 6);
        print(out, "alpha"  , 5);
        Table::print_init_header(out);
        print(out, "cap"    , 9);
        print(out, "n_full" , 9);
        print(out, "t_in"   , 8);
        print(out, "in_err" , 6);
        #ifdef MALLOC_COUNT
        print(out, "memory" , 7);
        #endif
        out << std::endl;
    }

    inline void print_timing(std::ostream& out,
                             size_t i,
                             size_t block,
                             double alpha,
                             size_t cap  ,
                             size_t n ,
                             double d_step,
                             size_t in_err,
                             Table& table)
    {
        print(out, i     , 4);
        print(out, block , 6);
        print(out, alpha , 5);
        table.print_init_data(out);
        print(out, cap   , 9);
        print(out, n     , 9);
        print(out, d_step, 8);
        print(out, in_err, 6);
        #ifdef MALLOC_COUNT
        print(out, double(malloc_count_current()-n*8)/double(8*2*block_size), 7);
        #endif
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t cap, size_t steps,
                   double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n; ++i)
        {
            keys[i] = dis(re);
        }

        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".reg",
                                      std::ofstream::out | std::ofstream::app);

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            Table table(cap, alpha, steps);

            auto in_errors = 0ull;

            for (size_t j = 0; j < n/block_size; ++j)
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                for (size_t ti = j*block_size; ti < (j+1)*block_size; ++ti)
                {
                    if (!table.insert(keys[ti], ti).second) ++in_errors;
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                double d_step = std::chrono::duration_cast<std::chrono::microseconds> (t2 - t1).count()/1000.;

                print_timing(*file, i, j, alpha, cap, n, d_step,
                         in_errors, table);
            }

        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    pin_to_core(0);
    CommandLine c(argn, argc);
    size_t      it    = c.intArg("-it"   , 5);
    size_t      n     = c.intArg("-n"    , 2000000);
    size_t      cap   = c.intArg("-cap"  , n);
    size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "");
    double      alpha = c.doubleArg("-alpha", 1.1);
    double      load  = c.doubleArg("-load" , 2.0);
    double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,no_hist_count> (c, it, n, cap, steps, alpha, name);
}
