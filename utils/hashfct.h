#ifndef HASHFCT_H
#define HASHFCT_H

struct crc_hasher
{
    inline size_t operator()(const long long& k) const
    {
        return size_t(__builtin_ia32_crc32di(1329235987123598723ull, k)
                   | (__builtin_ia32_crc32di(1383568923875084501ull, k) << 32));
    }
};


#endif // HASHFCT_H
