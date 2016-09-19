#include "include/spacegrow.h"
#include "utils/hashfct.h"

#include <random>
#include <iostream>




int main(int, char**)
{

    constexpr size_t n = 1000000;
    constexpr size_t range = (1ull<<63) -1;

    size_t* keys = new size_t[n];

    std::uniform_int_distribution<uint64_t> dis(1,range);
    std::mt19937_64 re;
    
    for (size_t i = 0; i < n; ++i)
    {
        keys[i] = dis(re);
    }
    
    std::cout << "randomness generated" << std::endl;

    auto errors = 0;
    bool first  = true;
    
    SpaceGrow<size_t, size_t, murmur_hasher> table(0, 1.1, 1007);
    
    std::cout << "table generated"      << std::endl;
    
    for (size_t i = 0; i < n; ++i)
    {
        if (!table.insert(keys[i], i))
        {
            ++errors;
            
            if (first)
            {
                first = false;
                std::cout << "first error at i=" << i << std::endl;
            }
        }
    }

    std::cout << "inserted elements encountered " << errors << " errors" << std::endl;

    auto count = 0;
    errors = 0;
    for (size_t i = 0; i < n; ++i)
    {
        auto e = table.find(keys[i]);
        if (e.first && (e.second == i)) count++;
        else if (e.first) std::cout << "wat" << std::endl;
        else errors++;
    }

    std::cout << "count: " << count << "  errors: " << errors << std::endl;

    table.printHist();

    delete[] keys;
    
    

    return 0;
}
