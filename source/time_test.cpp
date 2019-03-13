#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

#include "selection.h"
#include "utils/default_hash.h"
#include "utils/commandline.h"
#include "utils/thread_basics.h"

#ifdef MALLOC_COUNT
  #include "malloc_count.h"
  static constexpr bool malloc_mode = true;
  size_t get_malloc () { return malloc_count_current(); }
#else
  static constexpr bool malloc_mode = false;
  constexpr size_t get_malloc () { return 0; }
#endif

#ifdef RSS_COUNT
  #include <stdio.h>
  static constexpr bool rss_mode = true;
  size_t get_rss()
  {
      long rss = 0L;
      FILE* fp = NULL;
      if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
          return (size_t)0L;		/* Can't open? */
      if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
      {
          fclose( fp );
          return (size_t)0L;		/* Can't read? */
      }
      fclose( fp );
      return (size_t)rss;
  }
#else
  static constexpr bool rss_mode = false;
  constexpr size_t get_rss() { return 0; }
#endif


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
        if constexpr (malloc_mode) print(out, "memory" , 7);
        if constexpr (rss_mode)    print(out, "rss"    , 7);
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
                             size_t rss,
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
        if constexpr (malloc_mode)
            print(out, double(get_malloc())/double(8*2*n)-1., 7);
        if constexpr (rss_mode)
            print(out, rss, 7);
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

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            size_t start_rss = get_rss();

            Table table(cap, alpha, steps);

            auto in_errors = 0ull;
            auto fin_errors = 0ull;

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < n0 && in_errors < 100; ++i)
            {
                if (!table.insert(keys[i], i).second) return ++in_errors;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            for (size_t i = n0; i < n && in_errors < 100; ++i)
            {
                if (!table.insert(keys[i], i).second) ++in_errors;
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            size_t final_rss = get_rss() - start_rss;

            auto t2i = std::chrono::high_resolution_clock::now();
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

            double d_in0 = std::chrono::duration_cast<std::chrono::microseconds>
                (t1 - t0).count()/1000.;
            double d_in1 = std::chrono::duration_cast<std::chrono::microseconds>
                (t2 - t1).count()/1000.;
            double d_fn0 = std::chrono::duration_cast<std::chrono::microseconds>
                (t3 - t2i).count()/1000.;
            double d_fn1 = std::chrono::duration_cast<std::chrono::microseconds>
                (t4 - t3).count()/1000.;

            print_timing(*file, i, alpha, cap, n0, n,
                         d_in0, d_in1, d_fn0, d_fn1,
                         in_errors, fin_errors, final_rss, table);
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
    size_t      n0    = c.intArg("-pre"  , n/2);
    if (n0 > n) std::cout << "n0 (pre) has to be smaller than n! set n = n0 = " << n0 << std::endl;
    size_t      cap   = c.intArg("-cap"  , n);
    size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "");
    double      alpha = c.doubleArg("-alpha", 1.1);
    double      load  = c.doubleArg("-load" , 2.0);
    double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,false> (c, it, n, n0, cap, steps, alpha, name);
}
