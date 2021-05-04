#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

#include "selection.h"
#include "utils/default_hash.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/pin_thread.hpp"
#include "utils/output.hpp"
namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

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

template<class Config>
struct Test
{
    //using table_type = ProbIndependentBase<HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config> >;
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    int operator()(size_t it, size_t n,  size_t cap, size_t steps,
                   double alpha)
    {
        otm::out() << otm::width(4) << "# it"
                   << otm::width(8) << "alpha";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width(9) << "cap"
                   << otm::width(9) << "n_full"
                   << otm::width(10) << "t_in"
                   << otm::width(10) << "t_find+"
                   << otm::width(10) << "t_find-"
                   << otm::width(9) << "in_err"
                   << otm::width(9) << "fi_err";
        if constexpr (malloc_mode) otm::out() << otm::width(7) << "memory";
        if constexpr (rss_mode)    otm::out() << otm::width(7) << "rss";
        otm::out() << std::endl;

        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[2*n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < 2*n; ++i)
        {
            keys[i] = dis(re);
        }

        for (size_t i = 0; i < it; ++i)
        {
            size_t start_rss = get_rss();

            table_type table(cap, alpha, steps);

            auto in_errors = 0ull;
            auto fin_errors = 0ull;

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < n && in_errors < 100; ++i)
            {
                if (!table.insert(keys[i], i).second) ++in_errors;
            }
            auto t1 = std::chrono::high_resolution_clock::now();

            [[maybe_unused]] size_t final_rss = get_rss() - start_rss;

            auto t2 = std::chrono::high_resolution_clock::now();
            //const table_type& ctable = table;
            for (size_t i = 0; i < n; ++i)
            {
                auto e = table.find(keys[i]);
                if ( (e == table.end()) || ((*e).second != i))
                    //if (ctable.at(keys[i]) != i)
                {  fin_errors++;  }
            }
            auto t3 = std::chrono::high_resolution_clock::now();
            for (size_t i = n; i < 2*n; ++i)
            {
                auto e = table.find(keys[i]);
                if ( (e != table.end()) && (keys[(*e).second] != keys[i]))
                {  fin_errors++;  }
            }
            auto t4 = std::chrono::high_resolution_clock::now();

            // double d_in0 = std::chrono::duration_cast<std::chrono::microseconds>
            //     (t1 - t0).count()/1000.;
            double d_in  = std::chrono::duration_cast<std::chrono::microseconds>
                (t1 - t0).count()/1000.;
            double d_fn0 = std::chrono::duration_cast<std::chrono::microseconds>
                (t3 - t2).count()/1000.;
            double d_fn1 = std::chrono::duration_cast<std::chrono::microseconds>
                (t4 - t3).count()/1000.;

            otm::out() << otm::width(4) << i
                       << otm::width(8) << alpha;
            table.print_init_data(otm::out());
            otm::out() << otm::width(9) << cap
                       << otm::width(9) << n
                       << otm::width(10) << d_in
                       << otm::width(10) << d_fn0
                       << otm::width(10) << d_fn1
                       << otm::width(9) << in_errors
                       << otm::width(9) << fin_errors;
            if constexpr (malloc_mode)
                otm::out() << otm::width( 7)
                           << double(get_malloc())/double(8*2*n)-1.;
            if constexpr (rss_mode)
                otm::out() << otm::width( 7) << final_rss;
            otm::out() << std::endl;
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t      it    = c.int_arg("-it"   , 5);
    size_t      n     = c.int_arg("-n"    , 2000000);
    size_t      cap   = c.int_arg("-cap"  , n);
    size_t      steps = c.int_arg("-steps", 512);

    double      alpha = c.double_arg("-alpha", 1.1);
    double      load  = c.double_arg("-load" , 2.0);
    double      eps   = c.double_arg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".time";
        otm::out().set_file(name);
    }

    return Chooser::execute<Test,false> (c, it, n, cap, steps, alpha);
}
