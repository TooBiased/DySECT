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

#include <iostream>
#include <fstream>
#include <chrono>


// TWO PROBLEMS!! cutoff words at the end of the buffer
//                what happens if read reads only half a buffer
template<class Config>
struct test_type
{
    using table_type = HASHTYPE<size_t, size_t, utm::hash_tm::default_hash, Config>;
    static constexpr size_t bsize = 2*1024*1024;
    static constexpr size_t hseed = 13358259232739045019ull;

    int operator()(size_t it, size_t cap, double alpha, size_t steps,
                   std::string inf_name)
    {
        otm::out() << otm::width(4)  << "# it"
                   << otm::width(8)  << "alpha";
        table_type::print_init_header(otm::out());
        otm::out() << otm::width(9)  << "cap"
                   << otm::width(9)  << "time"
                   << otm::width(10) << "doubles"
                   << otm::width(9)  << "individ"
                   << otm::width(9)  << "err"
#ifdef MALLOC_COUNT
                   << otm::width(7)  << "mem"
                   << otm::width(7)  << "mx_mem"
#endif
                   << std::endl;


        for (size_t i = 0; i < it; ++i)
        {
            #ifdef MALLOC_COUNT
            malloc_count_reset_peak();
            #endif

            table_type table(cap, alpha, steps);
            std::ifstream inf (inf_name, std::ifstream::binary);
            char buffer[bsize];
            size_t contained  = 0;
            size_t individual = 0;
            size_t max        = 0;
            size_t errors     = 0;

            auto t0 = std::chrono::high_resolution_clock::now();
            while (inf.tellg() >= 0)
            {
                inf.read(buffer, bsize);

                char* j0 = buffer;
                for (char* j = buffer; j < buffer+inf.gcount(); ++j)
                {
                    if (*j == ' ' || *j == '\n')
                   {
                        if (j0 < j)
                        {
                            size_t key = XXH64(j0, j-j0, hseed);
                            auto e = table.insert(key, 1);
                            if (e.second) ++individual;
                            else
                            {
                                if (e.first != table.end())
                                {
                                    size_t a = ++((*e.first).second);
                                    max = (a < max) ? max : a;
                                    ++contained;
                                }
                                else ++errors;
                            }
                        }
                        j0 = j+1;
                    }
                }
                inf.seekg(j0-buffer-bsize, std::ios_base::cur);

            }
            auto t1 = std::chrono::high_resolution_clock::now();

            double ttest = std::chrono::duration_cast<std::chrono::microseconds>
                               (t1 - t0).count()/1000.;

            otm::out() << otm::width(4)  << i
                       << otm::width(8)  << alpha;
            table.print_init_data(otm::out());
            otm::out() << otm::width(9)  << cap
                       << otm::width(9)  << ttest
                       << otm::width(10) << contained
                       << otm::width(9)  << individual
                       << otm::width(9)  << errors
#ifdef MALLOC_COUNT
                       << otm::width(7)
                       << double(malloc_count_current())/double(8*2*n)
                       << otm::width(7)
                       << double(malloc_count_peak())/double(8*2*n)
#endif
                       << std::endl;

            inf.close();
        }

        return 0;
    }
};

int main(int argn, char** argc)
{
    utm::pin_to_core(0);
    utm::command_line_parser c(argn, argc);

    size_t      it    = c.int_arg("-it"   , 5);
    size_t      cap   = c.int_arg("-cap"  , 1000);
    size_t      steps = c.int_arg("-steps", 512);

    double      alpha = c.double_arg("-alpha", 1.1);
    double      load  = c.double_arg("-load" , 2.0);
    double      eps   = c.double_arg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    if (c.bool_arg("-out") || c.bool_arg("-file"))
    {
        std::string name = c.str_arg("-out", "");
        name = c.str_arg("-file", name) + ".crawl";
        otm::out().set_file(name);
    }

    if (c.bool_arg("-help") || c.bool_arg("-h"))
    {
        otm::out() << "This test measures access speed under a real world\n"
                   << "distribution, i.e., word count the file in -in param.\n"
                   << "\n"
                   << otm::color::blue << "1. " << otm::color::reset
                   << "hash each word, and insert the hash into the table\n"
                   << otm::color::blue << "2. " << otm::color::reset
                   << "if it was already present, increment its counter"
                   << std::endl;
        return 0;
    }

    const std::string inf   = c.str_arg("-in"   , "/home/maier/WorkEnv/tobias.org");

    return Chooser::execute<test_type, false>(c, it, cap, alpha, steps, inf);
}
