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
    //using table_type = ProbIndependentBase<HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config> >;
    using table_type = HASHTYPE<size_t, size_t, utils_tm::hash_tm::default_hash, Config>;

    int operator()(size_t it, size_t n, size_t pre, size_t cap, size_t pattern, size_t steps,
                   double alpha)
    {

        otm::out() << otm::width(4) <<  "# it"
                   << otm::width(8) <<  "alpha";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width(9) << "cap"
                   << otm::width(4) << "pat"
                   << otm::width(9) << "pre"
                   << otm::width(9) << "n"
                   << otm::width(8) << "t_pre"
                   << otm::width(8) << "t_mix"
                   << otm::width(9) << "in_err"
#ifdef MALLOC_COUNT
                   << otm::width(7) << "memory"
#endif
                   << std::endl;


        constexpr size_t range = (1ull<<63) -1;

        size_t  p_rep = n/10;
        size_t* ikeys = new size_t[pre+p_rep*pattern];
        size_t* fkeys = new size_t[p_rep*(10-pattern)];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < pre; ++i)
        {
            ikeys[i] = dis(re);
        }
        size_t pcount = 0;
        size_t icount = pre-1;
        size_t fcount = 0;
        for (size_t i = pre; i < p_rep*10+pre; ++i)
        {
            if (pcount < pattern)
            {
                ikeys[++icount] = dis(re);
            }
            else
            {
                std::uniform_int_distribution<uint64_t> find_dis(0,icount);
                size_t k = find_dis(re);
                fkeys[fcount++] = ikeys[k];
            }
            pcount  = (pcount < 9) ? pcount+1 : 0;
        }

        for (size_t i = 0; i < it; ++i)
        {
            table_type table(cap, alpha, steps);

            auto errors  = 0ull;
            auto ferrors = 0ull;

            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < pre; ++i)
            {
                if (! table.insert(ikeys[i], i).second) ++errors;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            pcount = 0;
            icount = pre;
            fcount = 0;
            for (size_t i = pre; i < p_rep*10+pre; ++i)
            {
                if (pcount < pattern)
                {
                    if (! table.insert(ikeys[icount++], i).second) ++errors;
                }
                else
                {
                    auto temp = table.find(fkeys[fcount++]);
                    if (temp == table.end()) ++ferrors;
                }
                pcount  = (pcount < 9) ? pcount+1 : 0;
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            double d_pre = std::chrono::duration_cast<std::chrono::microseconds> (t1 - t0).count()/1000.;
            double d_mix = std::chrono::duration_cast<std::chrono::microseconds> (t2 - t1).count()/1000.;

            otm::out() << otm::width(4) <<  i
                       << otm::width(8) <<  alpha ;
            table.print_init_data(otm::out());
            otm::out() << otm::width(9) <<  cap
                       << otm::width(4) <<  pattern
                       << otm::width(9) <<  pre
                       << otm::width(9) <<  n
                       << otm::width(8) <<  d_pre
                       << otm::width(8) <<  d_mix
                       << otm::width(9) <<  errors
#ifdef MALLOC_COUNT
                       << otm::width(7)
                       << double(malloc_count_current()-(n+pre)*8)/
                          (0.1*double(pattern*8*2*n)+double(8*2*pre))
#endif
                       << std::endl;

            if (ferrors && (!errors)) otm::out() << "Critical Error" << std::endl;
        }

        delete[] ikeys;
        delete[] fkeys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t      it    = c.int_arg("-it"   , 5);
    size_t      n     = c.int_arg("-n"    , 1000000);
    size_t      n0    = c.int_arg("-pre"  , n/4);
    size_t      cap   = c.int_arg("-cap"  , n0);
    size_t      pattern = c.int_arg("-pattern", 5);
    size_t      steps = c.int_arg("-steps", 512);

    double      alpha = c.double_arg("-alpha", 1.1);
    double      load  = c.double_arg("-load" , 2.0);
    double      eps   = c.double_arg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".mix";
        otm::out().set_file(name);
    }

    return Chooser::execute<Test,false> (c, it, n, n0, cap, pattern, steps, alpha);
}
