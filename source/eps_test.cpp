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

    constexpr static size_t block_size = 100000;

    inline void print_headline(std::ostream& out)
    {
        print(out, "# it"   , 4);
        Table::print_init_header(out);
        print(out, "n_step" , 9);
        print(out, "win"    , 5);
        print(out, "alpha"  , 6);
        print(out, "t_in"   , 8);
        print(out, "t_fip"  , 8);
        print(out, "t_fin"  , 8);
        print(out, "in_err" , 6);
        #ifdef MALLOC_COUNT
        print(out, "memory" , 7);
        #endif
        out << std::endl;
    }

    inline void print_timing(std::ostream& out,
                             size_t i,
                             size_t n,
                             size_t win,
                             double alpha,
                             double tin,
                             double tfip,
                             double tfin,
                             size_t in_err,
                             size_t fi_err,
                             size_t keysize,
                             Table& table)
    {
        print(out, i     , 4);
        table.print_init_data(out);
        print(out, n     , 9);
        print(out, win   , 5);
        print(out, alpha , 6);
        print(out, tin   , 8);
        print(out, tfip  , 8);
        print(out, tfin  , 8);
        print(out, in_err, 6);
        print(out, fi_err, 6);
        #ifdef MALLOC_COUNT
        print(out, malloc_count_current() - keysize*8, 7);
        #endif
        out << std::endl;
    }

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
                   double alpha, double alstp, size_t win, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t keysize = double(n)*alpha;
        size_t* inkeys = new size_t[keysize];
        size_t* fikeys = new size_t[win<<1];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".eps",
                                      std::ofstream::out | std::ofstream::app);

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            // prepare execution randomize each iteration seperately
            // to take variance into account
            for (size_t i = 0; i < keysize; ++i)
            {
                inkeys[i] = dis(re);
            }
            Table table(keysize, 1.0, steps);

            size_t in_errors = 0ull;
            size_t fi_errors = 0ull;

            size_t j       = 0;

            size_t nxt_blk = n;
            double ta      = alpha;
            prepare_find_keys(fikeys, inkeys, win, n);

            while (nxt_blk + win < keysize)
            {
                for (; j < nxt_blk; ++j)
                {
                    if (!table.insert(inkeys[j], j)) ++in_errors;
                }

                auto t1 = std::chrono::high_resolution_clock::now();
                for ( ; j < nxt_blk+win; ++j)
                {
                    if (!table.insert(inkeys[j], j)) ++in_errors;
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                for (size_t k = 0; k < win; ++k)
                {
                    auto temp = table.find(fikeys[k]);
                    if (!temp.first ) ++fi_errors;
                }
                auto t3 = std::chrono::high_resolution_clock::now();
                for (size_t k = win; k < win<<1; ++k)
                {
                    auto temp = table.find(fikeys[k]);
                    if ( temp.first ) ++fi_errors;
                }
                auto t4 = std::chrono::high_resolution_clock::now();
                double d_in  = std::chrono::duration_cast<std::chrono::nanoseconds> (t2 - t1).count()/double(win);
                double d_fip = std::chrono::duration_cast<std::chrono::nanoseconds> (t3 - t2).count()/double(win);
                double d_fin = std::chrono::duration_cast<std::chrono::nanoseconds> (t4 - t3).count()/double(win);

                if (in_errors > 100) break;
                print_timing(*file, i, nxt_blk, win, ta, d_in, d_fip, d_fin,
                             in_errors, fi_errors, keysize, table);

                prepare_find_keys(fikeys, inkeys, win, n);

                ta -= alstp;
                nxt_blk = double(keysize)/ta;
            }

        }

        delete[] inkeys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      it    = c.intArg("-it"      , 5);
    const size_t      n     = c.intArg("-n"       , 2000000);
    const size_t      steps = c.intArg("-steps"   , 512);
    const std::string name  = c.strArg("-out"     , "");
    const double      alpha = c.doubleArg("-alpha", 1.2);
    const double      alstp = c.doubleArg("-alsteps", 0.005);
    const size_t      win   = c.intArg("-win"     , 1000);

    return Chooser::execute<Test,no_hist_count>
        (c, it, n, steps, alpha, alstp, win, name);
}
