#include "selection.h"
#include "utils/hashfct.h"
#include "utils/commandline.h"
#include "utils/thread_basics.h"

#include <random>
#include <iostream>
#include <fstream>

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

    void print_hist(Table& table, std::string name)
    {
        std::ofstream out(name, std::ofstream::out);

        auto& hcount = table.hcounter;

        print(out, "# steps", 7);
        print(out, "nFitted", 8);
        out << std::endl;

        for (size_t i = 0; i < hcount.steps; ++i)
        {
            print(out, i, 7);
            print(out, hcount.hist[i], 8);
            out << std::endl;
        }
        out.close();
    }

    void print_dist(Table& table, std::string name)
    {
        std::ofstream out(name, std::ofstream::out);

        print (out, "# tab", 5);
        size_t gHist[table.bs+1];
        for (size_t i = 0; i <= table.bs; ++i)
        {
            gHist[i] = 0;
            print (out, i, 6);
        }

        print (out, "n"    , 8);
        print (out, "cap"  , 8);
        out << std::endl;

        for (size_t tl = 0; tl < table.tl; ++tl)
        {
            size_t lHist[table.bs+1];
            for (size_t i = 0; i <= table.bs; ++i) lHist[i] = 0;

            auto ltab = table.getTable(tl);

            for (size_t j = 0; j < ltab.first; ++j)
            {
                auto a = ltab.second[j].probe(0);
                if (a >= 0) ++lHist[a];
            }

            size_t n = 0;
            print (out, tl, 5);
            for (size_t i = 0; i <= table.bs; ++i)
            {
                print (out, lHist[i], 6);
                n += lHist[i] * (table.bs - i);
                gHist[i] += lHist[i];
            }
            print (out, n, 8);
            print (out, ltab.first*table.bs, 8);
            out << std::endl;
        }

        size_t n = 0;
        print (out, "#all", 5);
        for (size_t i = 0; i <= table.bs; ++i)
        {
            print (out, gHist[i], 6);
            n += gHist[i] * (table.bs - i);
        }
        print (out, n, 8);
        //print (out, table.curGrowAmount * (table.tl+table.curGrowTable), 8);
        out << std::endl;

        out.close();
    }

    int operator()(size_t n, size_t pre, size_t cap, size_t steps, double alpha, std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n+pre];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n + pre; ++i)
        {
            keys[i] = dis(re);
        }

        Table table(cap, alpha, steps);

        for (size_t i = 0; i < pre; ++i)
        {
            table.insert(keys[i], i);
        }

        table.clearHist();

        for (size_t i = pre; i < pre+n; ++i)
        {
            table.insert(keys[i], i);
        }

        print_dist(table, name + ".dist");
        print_hist(table, name + ".hist");

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    pin_to_core(0);
    CommandLine c(argn, argc);
    size_t      n     = c.intArg("-n"    , 1000000);
    size_t      pre   = c.intArg("-pre"  , 0);
    size_t      cap   = c.intArg("-cap"  , 1000000);
    size_t      steps = c.intArg("-steps", 512);
    std::string name  = c.strArg("-out"  , "temp");
    double      alpha = c.doubleArg("-alpha", 1.1);
    double      load  = c.doubleArg("-load" , 2.0);
    double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test, hist_count>(c, n, pre, cap, steps, alpha, name);
}
