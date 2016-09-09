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
    
    SpaceGrow<size_t, size_t, crc_hasher> table(n);
    
    std::cout << "table generated"      << std::endl;
    for (size_t i = 0; i < n; ++i)
    {
        std::cout << "i: " << i << std::endl;
        if (!table.insert(keys[i], i))
        {
            std::cout << "!!!" << std::endl;
            ++errors;
            
            if (first)
            {
                first = false;
                std::cout << "first error at i=" << i << std::endl;
            }
        }
        std::cout << "ei: " << i << std::endl;
    }

    delete[] keys;
    
    

    return 0;
}
