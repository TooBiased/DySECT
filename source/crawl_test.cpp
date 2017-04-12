#include "selection.h"
#include "utils/hashfct.h"
#include "utils/commandline.h"

#ifdef MALLOC_COUNT
#include "malloc_count.h"
#endif

#include <iostream>
#include <fstream>
#include <chrono>

template <class T>
inline void print(std::ostream& out, const T& t, size_t w)
{
    out.width(w);
    out << t << " " << std::flush;
}


// TWO PROBLEMS!! cutoff words at the end of the buffer
//                what happens if read reads only half a buffer
template<class Config>
struct Test
{
    using Table    = HASHTYPE<size_t, size_t, HASHFCT, Config>;
    static constexpr size_t bsize = 2*1024*1024;
    static constexpr size_t hseed = 13358259232739045019ull;


    void print_headline(std::ostream& out)
    {
        print(out, "# it"   , 4);
        print(out, "alpha"  , 6);
        Table::print_init_header(out);
        print(out, "cap"    , 9);
        print(out, "time"   , 9);
        print(out, "doubles", 10);
        print(out, "individ", 9);
        print(out, "err" , 6);
        #ifdef MALLOC_COUNT
        print(out, "mem"    , 7);
        print(out, "mx_mem" , 7);
        #endif

        out << std::endl;
    }

    void print_timing(std::ostream& out,
                      size_t it,
                      double alpha,
                      size_t cap,
                      double ttest,
                      size_t contained,
                      size_t n,
                      size_t errors,
                      Table& table)
    {
        print(out, it       , 4);
        print(out, alpha    , 6);
        table.print_init_data(out);
        print(out, cap      , 9);
        print(out, ttest    , 9);
        print(out, contained, 10);
        print(out, n        , 9);
        print(out, errors   , 6);
        #ifdef MALLOC_COUNT
        print(out, double(malloc_count_current())/double(8*2*n), 7);
        print(out, double(malloc_count_peak()   )/double(8*2*n), 7);
        #endif
        out << std::endl;
    }

    int operator()(size_t it, size_t cap, double alpha, size_t steps,
                   std::string inf_name, std::string outf_name)
    {
        std::ostream *outf;
        if (outf_name == "") outf = &(std::cout);
        else outf = new std::ofstream(outf_name + ".crawl",
                                      std::ofstream::out | std::ofstream::app);
        print_headline(*outf);

        for (size_t i = 0; i < it; ++i)
        {
            #ifdef MALLOC_COUNT
            malloc_count_reset_peak();
            #endif

            Table table(cap, alpha, steps);
            std::ifstream inf (inf_name, std::ifstream::binary);
            char buffer[bsize];
            size_t contained  = 0;
            size_t individual = 0;
            size_t max        = 0;
            size_t errors     = 0;

            if (! inf.is_open() )
            { std::cout << "file not found" << std::endl; return 5; }

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

            print_timing(*outf, i, alpha, cap, ttest, contained, individual, errors, table);

            inf.close();
        }

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      it    = c.intArg("-it"   , 5);
    const size_t      cap   = c.intArg("-cap"  , 1000);
    const size_t      steps = c.intArg("-steps", 512);
    const double      alpha = c.doubleArg("-alpha", 1.1);
    const double      load  = c.doubleArg("-load" , 2.0);
    const double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    const std::string inf   = c.strArg("-in"   , "/home/maier/WorkEnv/tobias.org");
    const std::string outf  = c.strArg("-out"  , "");

    return Chooser::execute<Test, hist_count>(c, it, cap, alpha, steps, inf, outf);
}
