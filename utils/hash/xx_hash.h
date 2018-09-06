#pragma once

// this define ensures, that xxhash is inlined/does not create new compile unit
#define XXH_PRIVATE_API
#include "xxhash.h"


namespace dysect {
namespace hash {

struct xx_hash
{
    xx_hash(size_t s = 13358259232739045019ull) : seed(s) { }

    static constexpr size_t significant_digits = 64;
    //const
    size_t seed;

    inline uint64_t operator()(const uint64_t k) const
    {
        auto local = k;
        return XXH64 (&local, 8, seed);
    }
};

}
}
