#pragma once

// We include the cpp to avoid generating another compile unit
#include "MurmurHash3.cpp"


namespace dysect {
namespace hash {

struct murmur3_hash
{
    murmur3_hash(size_t s = 1203989050u) : seed(s) { }

    static constexpr size_t significant_digits = 64;
    //const
    uint seed;

    inline uint64_t operator()(const uint64_t k) const
    {
        uint64_t local = k;
        uint64_t target[2];

        MurmurHash3_x64_128 (&local, 8, seed, target);

        return target[0];
    }
};

}
}
