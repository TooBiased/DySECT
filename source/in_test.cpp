#include "selection.h"
#include "utils/default_hash.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/pin_thread.hpp"
#include "utils/output.hpp"

#include <random>
#include <iostream>
#include <fstream>

namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

template<class Config>
struct test_type
{
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    void print_hist(table_type& table)
    {
        std::ofstream out(name, std::ofstream::out);

        auto& hcount = table.hcounter;

        otm::out() << otm::width(7) << "# steps"
                   << otm::width(8) << "nFitted"
                   << std::endl;

        for (size_t i = 0; i < hcount.steps; ++i)
        {
            otm::out() << otm::width(7) << i
                       << otm::width(8) << hcount.hist[i]
                       << std::endl;
        }
    }

    void print_dist(table_type& table)
    {
        otm::out() << otm::width(5) << "# tab";
        size_t gHist[table.bs+1];
        for (size_t i = 0; i <= table.bs; ++i)
        {
            gHist[i] = 0;
            otm::out() << otm::width(6) << i;
        }

        otm::out() << otm::width(8) << "n"
                   << otm::width(8) << "cap"
                   << std::endl;

        for (size_t tl = 0; tl < table.tl; ++tl)
        {
            size_t lHist[table.bs+1];
            for (size_t i = 0; i <= table.bs; ++i) lHist[i] = 0;

            auto ltab = table.gettable_type(tl);

            for (size_t j = 0; j < ltab.first; ++j)
            {
                auto a = ltab.second[j].probe(0);
                if (a >= 0) ++lHist[a];
            }

            size_t n = 0;
            otm::out() << otm::width(5) << tl;
            for (size_t i = 0; i <= table.bs; ++i)
            {
                otm::out() << otm::width(6) << lHist[i];
                n += lHist[i] * (table.bs - i);
                gHist[i] += lHist[i];
            }
            otm::out() << otm::width(8) << n
                       << otm::width(8) << ltab.first*table.bs
                       << std::endl;
        }

        size_t n = 0;
        otm::out() << otm::width(5) << "#all";
        for (size_t i = 0; i <= table.bs; ++i)
        {
            otm::out() <<  otm::width(6) << gHist[i];
            n += gHist[i] * (table.bs - i);
        }
        otm::out() <<  otm::width(8) << n
                   << std::endl;
    }

    int operator()(size_t n, size_t pre, size_t cap, size_t steps, double alpha)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n+pre];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n + pre; ++i)
        {
            keys[i] = dis(re);
        }

        table_type table(cap, alpha, steps);

        for (size_t i = 0; i < pre; ++i)
        {
            table.insert(keys[i], i);
        }

        table.clearHist();

        for (size_t i = pre; i < pre+n; ++i)
        {
            table.insert(keys[i], i);
        }

        otm::out().set_file(name + ".dist");
        print_dist(table);
        otm::out().set_file(name + ".hist");
        print_hist(table);

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t      n     = c.int_arg("-n"    , 1000000);
    size_t      pre   = c.int_arg("-pre"  , 0);
    size_t      cap   = c.int_arg("-cap"  , 1000000);
    size_t      steps = c.int_arg("-steps", 512);

    double      alpha = c.double_arg("-alpha", 1.1);
    double      load  = c.double_arg("-load" , 2.0);
    double      eps   = c.double_arg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    std::string name  = c.str_arg("-out"  , "temp");
    name              = c.str_arg("-file" , name);

    return Chooser::execute<test_type, hist_count>(c, n, pre, cap, steps, alpha, name);
}
