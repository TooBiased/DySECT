#include "selection.h"

#include "utils/default_hash.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/pin_thread.hpp"
#include "utils/output.hpp"
namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

#ifdef MALLOC_COUNT
#include "malloc_count.h"
#endif

#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

template<class Config>
struct test_type
{
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    int operator()(size_t it, size_t n, size_t win,  size_t cap, size_t steps,
                   double alpha)
    {
        otm::out() << otm::width(4) << "# it"
                   << otm::width(8) << "alpha";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width(9) << "cap"
                   << otm::width(9) << "n_win"
                   << otm::width(9) << "n_mix"
                   << otm::width(8) << "t_win"
                   << otm::width(8) << "t_mix"
                   << otm::width(8) << "t_eval"
                   << otm::width(9) << "op_err"
                   << otm::width(9) << "ev_err"
        #ifdef MALLOC_COUNT
                   << otm::width(7) << "mem"
        #endif
                   << std::endl;


        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n+win];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n+win; ++i)
        {
            keys[i] = dis(re);
        }

        for (size_t i = 0; i < it; ++i)
        {
            auto in_errors   = 0;
            auto del_errors  = 0;
            auto fin_errors  = 0;

            table_type table(cap, alpha, steps);

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0  ; i < win  ; ++i)
            {
                if (! table.insert(keys[i], i).second) ++in_errors;
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            for (size_t i = win; i < win+n; ++i)
            {
                if (!table.insert(keys[i], i).second)  ++in_errors;
                if (!table.erase(keys[i-win])) ++del_errors;
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0  ; i < n; ++i)
            {
                auto e = table.find(keys[i]);
                if (e != table.end())  ++fin_errors;;
            }
            for (size_t i = n; i < win+n; ++i)
            {
                auto e = table.find(keys[i]);
                if (e == table.end()) ++fin_errors;
            }

            for (size_t i = n; i < win+n; ++i)
            {
                if (!table.erase(keys[i])) ++del_errors;
            }

            auto t3 = std::chrono::high_resolution_clock::now();

            double d_win  = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1000.;
            double d_mix  = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()/1000.;
            double d_eval = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count()/1000.;

            otm::out() << otm::width(4) << i
                       << otm::width(8) << alpha;
            table.print_init_data(otm::out());
            otm::out() << otm::width(9) << cap
                       << otm::width(9) << win
                       << otm::width(9) << n
                       << otm::width(8) << d_win
                       << otm::width(8) << d_mix
                       << otm::width(8) << d_eval
                       << otm::width(7) << del_errors+in_errors
                       << otm::width(7) << fin_errors
#ifdef MALLOC_COUNT
                       << otm::width(7)
                       << double(malloc_count_current()-8*(n+win))/double(8*2*win)
#endif
                       << std::endl;
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
    size_t      n     = c.int_arg("-n"    , 1000000);
    size_t      win   = c.int_arg("-pre"  , n/2);
    size_t      cap   = c.int_arg("-cap"  , win);
    size_t      steps = c.int_arg("-steps", 512);

    double      alpha = c.double_arg("-alpha", 1.1);
    double      load  = c.double_arg("-load" , 2.0);
    double      eps   = c.double_arg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".del";
        otm::out().set_file(name);
    }

    if (c.bool_arg("-help") || c.bool_arg("-h"))
    {
        otm::out() << "This test measures access speed of deletions \n"
                   << "\n"
                   << otm::color::blue << "1. " << otm::color::reset
                   << "insert -pre random elements into a table with cap slots\n"
                   << otm::color::blue << "2. " << otm::color::reset
                   << "repeat -n insert+erase cycles\n"
                   << otm::color::blue << "3. " << otm::color::reset
                   << "evaluate all previous moves"
                   << std::endl;
        return 0;
    }

    return Chooser::execute<test_type,false> (c, it, n, win, cap, steps, alpha);
}
