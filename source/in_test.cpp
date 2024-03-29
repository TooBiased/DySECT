#include "selection.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/default_hash.hpp"
#include "utils/output.hpp"
#include "utils/pin_thread.hpp"

#include <fstream>
#include <iostream>
#include <random>

namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

static constexpr size_t TABLE_BS = 4;
static constexpr size_t TABLE_TL = 1;

class history_histogram
{
  public:
    history_histogram(size_t s) : steps(s), hist(new size_t[s])
    {
        for (size_t i = 0; i < s; ++i) { hist[i] = 0; }
    }

    void add(size_t i)
    {
        auto ind = (i < steps) ? i : steps - 1;
        ++hist[ind];
    }

    void clear()
    {
        for (size_t i = 0; i < steps; ++i) { hist[i] = 0; }
    }

    const size_t              steps;
    std::unique_ptr<size_t[]> hist;
};


template <class Config>
struct test_type
{
    using table_type =
        HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    void print_hist(table_type& table)
    {
        auto& hcount = table.get_history();

        otm::out() << otm::width(7) << "# steps" << otm::width(8) << "nFitted"
                   << std::endl;

        for (size_t i = 0; i < hcount.steps; ++i)
        {
            otm::out() << otm::width(7) << i << otm::width(8) << hcount.hist[i]
                       << std::endl;
        }
    }

    void print_dist(table_type& table)
    {
        otm::out() << otm::width(5) << "# tab";
        size_t gHist[TABLE_BS + 1];
        for (size_t i = 0; i <= TABLE_BS; ++i)
        {
            gHist[i] = 0;
            otm::out() << otm::width(6) << i;
        }

        otm::out() << otm::width(8) << "n" << otm::width(8) << "cap"
                   << std::endl;

        for (size_t tl = 0; tl < TABLE_TL; ++tl)
        {
            size_t lHist[TABLE_BS + 1];
            for (size_t i = 0; i <= TABLE_BS; ++i) lHist[i] = 0;

            auto ltab = table.getTable(tl);

            for (size_t j = 0; j < ltab.first; ++j)
            {
                auto a = ltab.second[j].probe(0);
                if (a >= 0) ++lHist[a];
            }

            size_t n = 0;
            otm::out() << otm::width(5) << tl;
            for (size_t i = 0; i <= TABLE_BS; ++i)
            {
                otm::out() << otm::width(6) << lHist[i];
                n += lHist[i] * (TABLE_BS - i);
                gHist[i] += lHist[i];
            }
            otm::out() << otm::width(8) << n << otm::width(8)
                       << ltab.first * TABLE_BS << std::endl;
        }

        size_t n = 0;
        otm::out() << otm::width(5) << "#all";
        for (size_t i = 0; i <= TABLE_BS; ++i)
        {
            otm::out() << otm::width(6) << gHist[i];
            n += gHist[i] * (TABLE_BS - i);
        }
        otm::out() << otm::width(8) << n << std::endl;
    }

    int operator()(size_t             n,
                   size_t             pre,
                   size_t             cap,
                   size_t             steps,
                   double             alpha,
                   const std::string& name)
    {
        constexpr size_t range = (1ull << 63) - 1;

        size_t* keys = new size_t[n + pre];

        std::uniform_int_distribution<uint64_t> dis(1, range);
        std::mt19937_64                         re;

        for (size_t i = 0; i < n + pre; ++i) { keys[i] = dis(re); }

        table_type table(cap, alpha, steps);

        for (size_t i = 0; i < pre; ++i) { table.insert(keys[i], i); }

        table.clear_history();

        for (size_t i = pre; i < pre + n; ++i) { table.insert(keys[i], i); }

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

    size_t n     = c.int_arg("-n", 1000000);
    size_t pre   = c.int_arg("-pre", 0);
    size_t cap   = c.int_arg("-cap", 1000000);
    size_t steps = c.int_arg("-steps", 512);

    double alpha = c.double_arg("-alpha", 1.1);
    double load  = c.double_arg("-load", 2.0);
    double eps   = c.double_arg("-eps", 1.0 - load);
    if (eps > 0.) alpha = 1. / (1. - eps);

    std::string name = c.str_arg("-out", "temp");
    name             = c.str_arg("-file", name);

    return Chooser::execute<test_type, history_histogram>(c, n, pre, cap, steps,
                                                          alpha, name);
}
