#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

#define EXTENDED
#include "selection.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/default_hash.hpp"
#include "utils/output.hpp"
#include "utils/pin_thread.hpp"
namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

class history_raw
{
  public:
    history_raw(size_t s = 0) { trend.reserve(s); }

    void add(size_t i) { trend.push_back(i); }
    void clear() { trend.clear(); }

    std::vector<size_t> trend;
};

template <class Config>
struct Test
{
    // using table_type = ProbIndependentBase<HASHTYPE<size_t, size_t,
    // dysect::hash::default_hash, Config> >;
    using table_type =
        HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    int operator()(size_t it, size_t n, size_t steps, size_t seed)
    {
        constexpr size_t range = (1ull << 63) - 1;

        size_t* keys = new size_t[n];

        for (size_t i = 0; i < it; ++i)
        {
            std::uniform_int_distribution<uint64_t> dis(1, range);
            std::mt19937_64                         re(seed * i);

            for (size_t j = 0; j < n; ++j) { keys[j] = dis(re); }

            table_type table(n, 1., steps);

            size_t ii = 0;

            auto t0 = std::chrono::high_resolution_clock::now();
            for (; ii < 0.995 * n; ++ii)
            {
                if (!table.insert(keys[ii], ii).second) break;
            }
            auto t1 = std::chrono::high_resolution_clock::now();

            otm::out() << "# seed:" << seed * i             //
                       << "  inserted:" << ii << "elements" //
                       << "  into:" << n << "slots"
                       << "  it took:"
                       << std::chrono::duration_cast<std::chrono::microseconds>(
                              t1 - t0)
                                  .count() /
                              1000.
                       << "ms" << std::endl;

            otm::out() << "#  element_nmbr    nmbr_of_slots    probed_slots    "
                          "inner_iteration"
                       << std::endl;

            auto& history = table.get_history();
            for (size_t j = history.trend.size() * 0.97; j < ii - 1; ++j)
            {
                otm::out() << j << "  " << n << "  " << history.trend[j]
                           << "  0" << std::endl;
            }

            for (size_t i = 0; i < it; ++i)
            {
                auto temptable = table.explicit_copy(n, 1., steps);

                for (size_t j = ii; j < n; ++j)
                {
                    if (!temptable.insert(dis(re), j).second) break;
                }
                auto& temphistory = temptable.get_history();
                for (size_t j = 0; j < temphistory.trend.size(); ++j)
                {
                    otm::out()
                        << ii + j << "  " << n << "  " << temphistory.trend[j]
                        << "  " << i << std::endl;
                }
            }
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t it    = c.int_arg("-it", 5);
    size_t n     = c.int_arg("-n", 2000000);
    size_t steps = c.int_arg("-steps", 512);
    size_t seed  = c.int_arg("-seed", 0xDEADBEEF);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name             = c.str_arg("-file", name) + ".probe";
        otm::out().set_file(name);
    }

    return Chooser::execute<Test, history_raw>(c, it, n, steps, seed);
}
