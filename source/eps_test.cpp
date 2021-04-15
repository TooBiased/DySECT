//#include "include/spacegrow.h"
#include "selection.h"

// #include "include/strategies/dstrat_bfs.h"
// #include "include/strategies/dstrat_rwalk.h"
// #include "include/strategies/dstrat_rwalk_cyclic.h"

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
    //using table_type = ProbIndependentBase<HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config> >;
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    constexpr static size_t block_size = 100000;

    void prepare_find_keys(size_t* fikeys, size_t* inkeys, size_t fis, size_t ins)
    {
        std::uniform_int_distribution<uint64_t> pdis(1,ins);
        std::mt19937_64 re(ins*2098740980990909098ull);
        for (size_t i = 0; i < fis; ++i)
        {
            fikeys[i] = inkeys[pdis(re)];
        }
        std::uniform_int_distribution<uint64_t> ndis(1,(1ull<<63) -1);
        for (size_t i = fis; i < fis<<1; ++i)
        {
            fikeys[i] = ndis(re);
        }
    }

    int operator()(size_t it, size_t n, size_t steps,
                   double eps, double epstp, size_t win)
    {
        otm::out() << otm::width(4) << "# it";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width(9) << "n_step"
                   << otm::width(5) << "win"
                   << otm::width(6) << "load"
                   << otm::width(8) << "t_in"
                   << otm::width(8) << "t_fip"
                   << otm::width(8) << "t_fin"
                   << otm::width(8) << "in_err"
                   << otm::width(8) << "fi_err"
#ifdef MALLOC_COUNT
                   << otm::width(7) << "memory"
#endif
                   << std::endl;

        constexpr size_t range = (1ull<<63) -1;

        size_t* inkeys = new size_t[n];
        size_t* fikeys = new size_t[win<<1];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < it; ++i)
        {

            // prepare execution randomize each iteration seperately
            // to take variance into account
            for (size_t i = 0; i < n; ++i)
            {
                inkeys[i] = dis(re);
            }
            table_type table(n, 1.0, steps);

            size_t in_errors = 0ull;
            size_t fi_errors = 0ull;

            size_t j       = 0;

            double tlf     = 1.-eps;
            size_t nxt_blk = n*tlf;

            prepare_find_keys(fikeys, inkeys, win, nxt_blk);

            while (nxt_blk + win < n)
            {
                for (; j < nxt_blk; ++j)
                {
                    if (!table.insert(inkeys[j], j).second) ++in_errors;
                }

                auto t1 = std::chrono::high_resolution_clock::now();
                for ( ; j < nxt_blk+win; ++j)
                {
                    if (!table.insert(inkeys[j], j).second) ++in_errors;
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                for (size_t k = 0; k < win; ++k)
                {
                    auto temp = table.find(fikeys[k]);
                    if (temp == table.end()) ++fi_errors;
                }
                auto t3 = std::chrono::high_resolution_clock::now();
                for (size_t k = win; k < win<<1; ++k)
                {
                    auto temp = table.find(fikeys[k]);
                    if (temp != table.end()) ++fi_errors;
                }
                auto t4 = std::chrono::high_resolution_clock::now();
                double d_in  = std::chrono::duration_cast<std::chrono::nanoseconds>
                    (t2 - t1).count()/double(win);
                double d_fip = std::chrono::duration_cast<std::chrono::nanoseconds>
                    (t3 - t2).count()/double(win);
                double d_fin = std::chrono::duration_cast<std::chrono::nanoseconds>
                    (t4 - t3).count()/double(win);

                if (in_errors > 100) break;

                otm::out() << otm::width(4) << i;
                table.print_init_data(otm::out());
                otm::out() << otm::width(9) << nxt_blk
                           << otm::width(5) << win
                           << otm::width(6) << tlf
                           << otm::width(8) << d_in
                           << otm::width(8) << d_fip
                           << otm::width(8) << d_fin
                           << otm::width(8) << in_errors
                           << otm::width(8) << fi_errors
#ifdef MALLOC_COUNT
                           << otm::width(7)
                           << malloc_count_current() - n*8 -win*16
#endif
                           << std::endl;

                tlf += epstp;
                nxt_blk = n*tlf;

                prepare_find_keys(fikeys, inkeys, win, nxt_blk);
            }

        }

        delete[] inkeys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t      it    = c.int_arg("-it"      , 5);
    size_t      n     = c.int_arg("-n"       , 2000000);
    size_t      steps = c.int_arg("-steps"   , 512);

    double      alpha = c.double_arg("-alpha", 0.1);
    double      load  = c.double_arg("-load" , 1.0-(alpha-1.0)/alpha);
    double      eps   = c.double_arg("-eps"  , 1.0-load);

    const double      epstp = c.double_arg("-epsteps", 0.005);
    const size_t      win   = c.int_arg("-win"     , 1000);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".eps";
        otm::out().set_file(name);
    }

    if (c.bool_arg("-help") || c.bool_arg("-h"))
    {
        otm::out() << "Test the performance depending on the fill degree\n"
                   << "stepwise test after blocks\n"
                   << "\n"
                   << otm::color::blue << "1. " << otm::color::reset
                   << "insert to the next fill degree test window\n"
                   << otm::color::blue << "2. " << otm::color::reset
                   << "test insert performance with -win new inserts\n"
                   << otm::color::blue << "3. " << otm::color::reset
                   << "query -win random inserted elements\n"
                   << otm::color::blue << "4. " << otm::color::reset
                   << "query -win not inserted elements\n"
                   << otm::color::blue << "repeat" << otm::color::reset
                   << std::endl;
        return 0;
    }

    if (eps < 0.) { std::cout << "please input eps" << std::endl; return 8;}

    return Chooser::execute<test_type,false>
        (c, it, n, steps, eps, epstp, win);
}
