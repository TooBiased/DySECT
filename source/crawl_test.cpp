#include "selection.h"
#include "utils/hashfct.h"
#include "utils/commandline.h"

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
        print(out, "it"     , 3);
        print(out, "alpha"  , 6);
        Table::print_init_header(out);
        print(out, "cap"    , 9);
        print(out, "time"   , 9);
        print(out, "doubles", 8);
        print(out, "individ", 8);
        out << std::endl;
    }

    void print_timing(std::ostream& out,
                      size_t it,
                      double alpha,
                      size_t cap,
                      double ttest,
                      size_t contained,
                      size_t n,
                      Table& table)
    {
        print(out, it       , 3);
        print(out, alpha    , 6);
        table.print_init_data(out);
        print(out, cap      , 9);
        print(out, ttest    , 9);
        print(out, contained, 8);
        print(out, n        , 8);
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
            Table table(cap, alpha, steps);
            std::ifstream inf (inf_name, std::ifstream::binary);
            char buffer[bsize];
            size_t contained = 0;

            if (! inf.is_open() )
            { std::cout << "file not found" << std::endl; return 5; }

            auto t0 = std::chrono::high_resolution_clock::now();
            while (inf.tellg() >= 0)
            {
                inf.read(buffer, bsize);
                //std::cout << inf.tellg() << " " << std::flush;
                //std::cout << inf.gcount() << std::endl;

                char* j0 = buffer;
                for (char* j = buffer; j < buffer+inf.gcount(); ++j)
                {
                    if (*j == ' ' || *j == '\n')
                    {
                        if (j0 < j)
                        {
                            size_t key = XXH64(j0, j-j0, hseed);
                            if (! table.insert(key, 1)) ++contained;
                            else
                            {
                                //outf->write(j0, j-j0);
                                //std::cout << " " << std::flush;
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

            print_timing(*outf, i, alpha, cap, ttest, contained, table.n, table);

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

    const std::string inf   = c.strArg("-in"   , "/home/maier/WorkEnv/tobias.org");
    const std::string outf  = c.strArg("-out"  , "");

    return Chooser::execute<Test, hist_count>(c, it, cap, alpha, steps, inf, outf);
}
