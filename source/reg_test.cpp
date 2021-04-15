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
struct Test
{
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    constexpr static size_t block_size = 100000;

    int operator()(size_t it, size_t n, size_t cap, size_t steps,
                   double alpha)
    {
        otm::out() << otm::width(4) << "# it"
                   << otm::width(6) << "block"
                   << otm::width(8) << "alpha";
        table_type::print_init_header(out);
        otm::out() << otm::width( 9) << "cap"
                   << otm::width( 9) << "n_full"
                   << otm::width( 8) << "t_in"
                   << otm::width( 6) << "in_err"
#ifdef MALLOC_COUNT
                   << otm::width( 7) << "memory"
#endif
                   << std::endl;


        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n; ++i)
        {
            keys[i] = dis(re);
        }

        for (size_t i = 0; i < it; ++i)
        {
            table_type table(cap, alpha, steps);

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

                otm::out() << otm::width(4) << i
                           << otm::width(6) << j
                           << otm::width(8) << alpha;
                table.print_init_data(out);
                otm::out() << otm::width(9) << cap
                           << otm::width(9) << n
                           << otm::width(8) << d_step
                           << otm::width(6) << in_errors
#ifdef MALLOC_COUNT
                           << otm::width(7) << double(malloc_count_current()-n*8)
                                               /double(8*2*block_size)
#endif
                           << std::endl;
            }
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utils_tm::pin_to_core(0);
    utils_tm::command_line_parser c(argn, argc);

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
        name = c.str_arg("-file", name) + ".reg";
        otm::out().set_file(name);
    }

    return Chooser::execute<Test,no_hist_count> (c, it, n, cap, steps, alpha);
}
