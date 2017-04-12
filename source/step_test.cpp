#include "selection.h"
#include "include/growing_cuckoo.h"
#include "include/simple_cuckoo.h"
#include "include/hom_2lvl_cuckoo.h"

#include "utils/hashfct.h"
#include "utils/commandline.h"

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
    using Grows  = GrowingCuckoo<size_t, size_t, HASHFCT, Config>;
    using Simple = SimpleCuckoo <size_t, size_t, HASHFCT, Config>;
    using Hom2Lvl= Hom2LvlCuckoo<size_t, size_t, HASHFCT, Config>;
    using MGrow  =

    template <class Table>
    void print_hist(Table& table, size_t lvl, std::string name)
    {
        std::ofstream out(name, std::ofstream::app);

        auto& hcount = table.hcounter;

        print(out, "# lvl", 5);
        print(out, "steps", 7);
        print(out, "nFitted", 8);
        out << std::endl;

        for (size_t i = 0; i < hcount.steps; ++i)
        {
            print(out, lvl, 5);
            print(out, i, 7);
            print(out, hcount.hist[i], 8);
            out << std::endl;
        }
        out.close();
    }

    template <class Table>
    void print_dist(Table& table, size_t lvl, std::string name)
    {
        std::ofstream out(name, std::ofstream::app);

        print (out, "# lvl", 5);
        print (out, "tab", 5);
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

            auto ctab = table.getTable(tl);

            for (size_t j = 0; j < ctab.first; ++j)
            {
                auto a = ctab.second[j].probe(0);
                if (a >= 0) ++lHist[a];
            }

            size_t n = 0;
            print (out, lvl, 5);
            print (out, tl, 5);
            for (size_t i = 0; i <= table.bs; ++i)
            {
                print (out, lHist[i], 6);
                n += lHist[i] * (table.bs - i);
                gHist[i] += lHist[i];
            }
            print (out, n, 8);
            print (out, ctab.first*table.bs, 8);
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

    void print_stuff(Grows& g_d, Grows& g_s,
                     Simple& s_s, Hom2Lvl& h_s,
                     size_t lvl, std::string name)
    {
        print_hist(g_d, lvl, name+"_grows_dynamic.hist");
        print_hist(g_s, lvl, name+"_grows_static.hist");
        print_hist(s_s, lvl, name+"_simple_static.hist");
        print_hist(h_s, lvl, name+"_hom2lvl_static.hist");

        print_dist(g_d, lvl, name+"_grows_dynamic.dist");
        print_dist(g_s, lvl, name+"_grows_static.dist");
        print_dist(s_s, lvl, name+"_simple_static.dist");
        print_dist(h_s, lvl, name+"_hom2lvl_static.dist");
    }

    void init_output_files(std::string name)
    {
        std::ofstream out0(name+"_grows_dynamic.hist", std::ofstream::out);
        out0.close();
        std::ofstream out1(name+"_grows_static.hist" , std::ofstream::out);
        out1.close();
        std::ofstream out2(name+"_simple_static.hist" , std::ofstream::out);
        out2.close();
        std::ofstream out3(name+"_hom2lvl_static.hist" , std::ofstream::out);
        out3.close();
        std::ofstream out4(name+"_grows_dynamic.dist", std::ofstream::out);
        out4.close();
        std::ofstream out5(name+"_grows_static.dist" , std::ofstream::out);
        out5.close();
        std::ofstream out6(name+"_simple_static.dist" , std::ofstream::out);
        out6.close();
        std::ofstream out7(name+"_hom2lvl_static.dist" , std::ofstream::out);
        out7.close();
    }

    template <class Table>
    void catchup(Table& table, size_t* keys, size_t n)
    {
        for (size_t i = 0; i < n; ++i)
        {
            table.insert(keys[i], i);
        }
    }



    int operator()(size_t n, size_t cap, size_t steps, double alpha, std::string name)
    {
        init_output_files(name);

        constexpr size_t range = (1ull<<63) -1;

        size_t* keys = new size_t[n];

        std::uniform_int_distribution<uint64_t> dis(1,range);
        std::mt19937_64 re;

        for (size_t i = 0; i < n; ++i)
        {
            keys[i] = dis(re);
        }

        Grows   g_d(cap             , alpha, steps);
        Grows   g_s(g_d.capacity, 1.0  , steps);
        Simple  s_s(g_d.capacity, 1.0  , steps);
        Hom2Lvl h_s(g_d.capacity, 1.0  , steps);

        auto curr_d_size = g_d.capacity;
        auto e_g_d = 0;
        auto e_g_s = 0;
        auto e_s_s = 0;
        auto e_h_s = 0;
        auto lvl = 0;

        for (size_t i = 0; i < n; ++i)
        {
            if (!g_d.insert(keys[i], i).second) { ++e_g_d; }
            if (!g_s.insert(keys[i], i).second) { ++e_g_s; }
            if (!s_s.insert(keys[i], i).second) { ++e_s_s; }
            if (!h_s.insert(keys[i], i).second) { ++e_h_s; }
            if ( g_d.capacity != curr_d_size )
            {
                print_stuff(g_d, g_s, s_s, h_s, lvl++, name);
                curr_d_size = g_d.capacity;

                g_s = std::move(Grows(curr_d_size, 1.0, steps));
                s_s = std::move(Simple(curr_d_size, 1.0, steps));
                h_s = std::move(Hom2Lvl(curr_d_size, 1.0, steps));

                std::cout << "NEW TABLES size: " << curr_d_size  << " "
                          << g_s.capacity << " " << s_s.capacity << " "
                          << h_s.capacity << std::endl;

                catchup(g_s, keys, i);
                catchup(s_s, keys, i);
                catchup(h_s, keys, i);

                g_d.clearHist();
                g_s.clearHist();
                s_s.clearHist();
                h_s.clearHist();
            }
        }

        std::cout << "inserted elements encountered " << e_g_d << " " << e_g_s
                  << " " << e_s_s << " " << e_h_s << " errors" << std::endl;

        delete[] keys;

        return 0;
    }
};

int main(int argn, char** argc)
{
    CommandLine c(argn, argc);
    const size_t      n     = c.intArg("-n"    , 1000000);
    const size_t      cap   = c.intArg("-cap"  , 10000);
    const size_t      steps = c.intArg("-steps", 512);
    const std::string name  = c.strArg("-out"  , "temp");
    const double      alpha = c.doubleArg("-alpha", 1.1);
    const double      load  = c.doubleArg("-load" , 2.0);
    const double      eps   = c.doubleArg("-eps"  , 1.0-load);
    if (eps > 0.) alpha = 1./(1.-eps);

    return Chooser::execute<Test,hist_count> (c, n, cap, steps, alpha, name);
}
