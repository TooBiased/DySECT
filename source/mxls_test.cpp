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

    int operator()(size_t it, size_t n, size_t steps)
    {
        otm::out() << otm::width( 4) << "# it"
                   << otm::width( 9) << "cap"
                   << otm::width( 9) << "n"
                   << otm::width( 6) << "steps"
                   << otm::width( 6) << "ngrow";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width( 9) << "g_elem"
                   << std::endl;

        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < it; ++i)
        {
            for (size_t i = 0; i < n; ++i)
            {
                keys[i] = dis(re);
            }

            table_type table(n, 1.0, steps);

            size_t ngrow = 0;
            for (size_t j = 0; j < n; ++j)
            {
                if (!table.insert(keys[j], j).second)
                {
                    otm::out() << otm::width(4) << i
                               << otm::width(9) << n
                               << otm::width(9) << n
                               << otm::width(6) << steps
                               << otm::width(6) << ngrow;
                    table.print_init_data(otm::out());
                    otm::out() << otm::width(9) << j
                               << std::endl;

                    ngrow++;
                    break;
                }
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

    size_t it     = c.int_arg("-it"   , 5);
    size_t n      = c.int_arg("-n"    , 2000000);
    size_t steps  = c.int_arg("-steps", 512);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".mxls";
        otm::out().set_file(name);
    }

    return Chooser::execute<Test,true> (c, it, n, steps);
}
