#pragma once

#include <cstddef>
#include <cstdint>


namespace dysect {
namespace hash {

struct crc_hash
{
    crc_hash(size_t seed = 12923598712359872066ull)
        : seed0(seed), seed1(seed*7467732452331123588ull)
    { }

    static const size_t significant_digits = 64;
    size_t seed0;
    size_t seed1;

    inline uint64_t operator()(const uint64_t& k) const
    {
        return uint64_t(   __builtin_ia32_crc32di(k, seed0)
                        | (__builtin_ia32_crc32di(k, seed1) << 32));
    }
};

}
}
