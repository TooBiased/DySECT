//#include "include/spacegrow.h"
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
    //using Table = ProbIndependentBase<HASHTYPE<size_t, size_t, HASHFCT, Config> >;
    using Table = HASHTYPE<size_t, size_t, HASHFCT, Config>;

    inline void print_headline(std::ostream& out)
    {
        print(out, "# it"   , 4);
        print(out, "cap"    , 9);
        print(out, "n"      , 9);
        print(out, "steps"  , 6);
        print(out, "ngrow"  , 6);
        Table::print_init_header(out);
        print(out, "g_elem" , 9);
        out << std::endl;
    }

    inline void print_timing(std::ostream& out,
                             size_t i,
                             size_t cap,
                             size_t n,
                             size_t steps,
                             size_t ngrow,
                             size_t j,
                             Table& table)
    {
        print(out, i     , 4);
        print(out, cap   , 9);
        print(out, n     , 9);
        print(out, steps , 6);
        print(out, ngrow , 6);
        table.print_init_data(out);
        print(out, j     , 9);
        out << std::endl;
    }

    int operator()(size_t it, size_t n, size_t cap, size_t steps,
                   std::string name)
    {
        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[2*n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        std::ostream* file;
        if (name == "") file = &(std::cout);
        else file = new std::ofstream(name + ".sing",
                                      std::ofstream::out | std::ofstream::app);

        print_headline(*file);

        for (size_t i = 0; i < it; ++i)
        {
            for (size_t i = 0; i < n; ++i)
            {
                keys[i] = dis(re);
            }

            Table table(cap, 1.0, steps);

            size_t ngrow = 0;
            for (size_t j = 0; j < n; ++j)
            {
                while (!table.insert(keys[j], j).second)
                {
                    print_timing(*file, i, cap, n, steps, ngrow++,  j, table);
                    table.explicit_grow();
                }
            }
        }

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    size_t it     = c.intArg("-it"   , 5);
    size_t n      = c.intArg("-n"    , 2000000);
    size_t cap    = c.intArg("-cap"  , 50000);
    size_t steps  = c.intArg("-steps", 512);

    const std::string name  = c.strArg("-out"  , "");

    return Chooser::execute<Test,true> (c, it, n, cap, steps, name);
}
