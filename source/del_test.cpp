#include "selection.h"

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
    using Table = HASHTYPE<size_t, size_t, HASHFCT, Config>;

    inline void print_headline(std::ostream& out)
    {
        print(out, "# it"  , 4);
        print(out, "alpha" , 5);
        Table::print_init_header(out);
        print(out, "cap"   , 9);
        print(out, "n_win" , 9);
        print(out, "n_mix" , 9);
        print(out, "t_win" , 8);
        print(out, "t_mix" , 8);
        print(out, "t_eval", 8);
        print(out, "op_err", 7);
        print(out, "ev_err", 7);
        #ifdef MALLOC_COUNT
        print(out, "mem"   , 7);
        #endif
        out << std::endl;
    }

    inline void print_timing(std::ostream& out, size_t i,
                             double alpha,
                             size_t cap,
                             size_t win,
                             size_t n,
                             double d_win,
                             double d_mix,
                             double d_eval,
                             size_t op_err,
                             size_t ev_err,
                             Table& table)
    {
        print(out, i     , 4);
        print(out, alpha , 5);
        table.print_init_data(out);
        print(out, cap   , 9);
        print(out, win   , 9);
        print(out, n     , 9);
        print(out, d_win , 8);
        print(out, d_mix , 8);
        print(out, d_eval, 8);
        print(out, op_err, 7);
        print(out, ev_err, 7);
        #ifdef MALLOC_COUNT
        print(out, double(malloc_count_current()-8*(n+win))/double(8*2*win), 7);
        #endif
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t win,  size_t cap, size_t steps,
                   double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n+win];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n+win; ++i)
        {
            keys[i] = dis(re);
        }

        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".del",
                                      std::ofstream::out | std::ofstream::app);

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            auto in_errors   = 0;
            auto del_errors  = 0;
            auto fin_errors  = 0;

            Table table(cap, alpha, steps);

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

            print_timing(*file, i, alpha, cap, win, n, d_win, d_mix, d_eval,
                         del_errors+in_errors, fin_errors, table);
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      it    = c.intArg("-it"   , 5);
    const size_t      n     = c.intArg("-n"    , 1000000);
    const size_t      win   = c.intArg("-pre"  , n/2);
    const size_t      cap   = c.intArg("-cap"  , win);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "");
    const double      alpha = c.doubleArg("-alpha", 1.1);
    const double      load  = c.doubleArg("-load" , 2.0);
    const double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,no_hist_count> (c, it, n, win, cap, steps, alpha, name);
}
