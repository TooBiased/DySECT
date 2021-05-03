
#include <random>
#include <limits>
#include <chrono>


#include "utils/default_hash.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/pin_thread.hpp"
#include "utils/output.hpp"
namespace utm = utils_tm;
namespace otm = utils_tm::out_tm;

#include "selection.h"

#ifdef MALLOC_COUNT
#include "malloc_count.h"
#endif


template<class Config>
struct test_type
{
    //using table = ProbIndependentBase<HASHTYPE<size_t, size_t, dysect::hash::default_hash, Config> >;
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;

    int operator()(size_t it, size_t n, size_t cap, size_t mdisp)
    {
        otm::out() << otm::width(3)  << "#it"
                   << otm::width(10) << "n"
                   << otm::width(10) << "cap"
                   << otm::width(5)  << "disp"
                   << otm::width(9)  << "ndisp"
                   << otm::width(10) << "time"
                   << std::endl;

        constexpr size_t range = std::numeric_limits<size_t>::max();

        size_t*              keys  = new size_t[n];
        std::vector<size_t>* disps = new std::vector<size_t>[mdisp+1];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t j = 0; j < it; ++j)
        {
            for (size_t i = 0; i < n; ++i)
            {
                keys[i] = dis(re);
            }
            for (size_t i = 0; i <= mdisp; ++i)
            {
                disps[i].clear();
            }

            table_type table(cap, 1.0, 1000);

            for (size_t i = 0; i < n; ++i)
            {
                if (!table.insert(keys[i], i).second)
                {
                    otm::out() << otm::color::red   << "missed insert key: "
                               << otm::color::reset << keys[i] << std::endl;
                }
            }

            for (size_t i = 0; i < n; ++i)
            {
                auto key  = keys[i];
                int  disp = std::min<int>(table.displacement(key), mdisp);

                if (disp < 0)
                {
                    otm::out() << otm::color::red   << "no displacement key: "
                               << otm::color::reset << key << std::endl;
                }

                disps[disp].push_back(key);
            }


            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < n; ++i)
            {
                if (table.find(keys[i]) == table.end())
                {
                    otm::out() << otm::color::red   << "unsuccessful query in cache clean"
                               << otm::color::reset << keys[i] << std::endl;
                }
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            double tfind = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()/1000.;

            for (size_t i = 0; i <= mdisp; ++i)
            {
                otm::out() << otm::width(3)  << j
                           << otm::width(10) << n
                           << otm::width(10) << cap
                           << otm::width(5)  << i
                           << otm::width(9)  << disps[i].size()
                           << std::flush;

                if (! disps[i].size())
                {
                    otm::out() << otm::width(10) << 0
                               << otm::width(10) << tfind
                               << std::endl;
                    continue;
                }

                // clean the cash
                for (size_t i = 0; i < n; ++i)
                {
                    if (table.find(keys[i]) == table.end())
                    {
                        otm::out() << otm::color::red   << "unsuccessful query in cache clean"
                                   << otm::color::reset << keys[i] << std::endl;
                    }
                }

                auto t0 = std::chrono::high_resolution_clock::now();
                for (auto& key : disps[i])
                {
                    if (table.find(key) == table.end())
                        otm::out() << otm::color::red << "unsuccessful query key: "
                                   << otm::color::reset << key << std::endl;
                }
                auto t1 = std::chrono::high_resolution_clock::now();
                double tdiff = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
                otm::out() << otm::width(10) << tdiff/disps[i].size()
                           << otm::width(10) << tfind << std::endl;
            }
        }


        delete[] keys;
        delete[] disps;

        return 0;
    }
};

int main(int argn, char** argc)
{
    // pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t it     = c.int_arg("-it" , 5);
    size_t n      = c.int_arg("-n"  , 700000);
    size_t cap    = c.int_arg("-cap", 1000000);
    size_t mdisp  = c.int_arg("-max_disp", 50);

    if (c.bool_arg("-file") || c.bool_arg("-out"))
    {
        std::string name = c.str_arg("-file", "");
        name = c.str_arg("-out", name) + ".displ";
        otm::out().set_file(name);
    }

    if (c.bool_arg("-help") || c.bool_arg("-h"))
    {
        otm::out() << "This test measures access speed dependant on the \n"
                   << "displacement of elements. For different hash tables\n"
                   << "\n"
                   << otm::color::blue << "1. " << otm::color::reset
                   << "insert n random elements into a table with cap slots\n"
                   << otm::color::blue << "2. " << otm::color::reset
                   << "for each element, measure its displacement\n"
                   << otm::color::blue << "3. " << otm::color::reset
                   << "iterate over displacements query all such elements\n"
                   << "   output their number, and the average time per query"
                   << std::endl;
        return 0;
    }

    return Chooser::execute<test_type,true> (c, it, n, cap, mdisp);
}
