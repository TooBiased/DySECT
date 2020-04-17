#pragma once

/*******************************************************************************
 * utils/default_hash.h
 *
 * In this file we chose an adequat default hash function, if nothing
 * else is defined XXHASH is chosen.
 *
 * If you have any problems with third party codes try defining MURMUR2.
 *
 * Part of Project growt - https://github.com/TooBiased/growt.git
 *
 * Copyright (C) 2015-2016 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#if (! (defined(CRC)     || \
        defined(MURMUR2) || \
        defined(MURMUR3) || \
        defined(XXHASH)) )
#define XXHASH
#endif // NO HASH DEFINED

#ifdef CRC
#include "hash/crc_hash.h"
#define HASHFCT dysect::hash::crc_hash
namespace dysect{
namespace hash {

    using default_hash = crc_hash;

}
}
#endif //CRC


#ifdef MURMUR2
#include "hash/murmur2_hash.h"
#define HASHFCT dysect::hash::murmur2_hash
namespace dysect{
namespace hash {

    using default_hash = murmur2_hash;

}
}
#endif // MURMUR2


#ifdef MURMUR3
#include "hash/murmur3_hash.h"
#define HASHFCT dysect::hash::murmur3_hash
namespace dysect{
namespace hash {

    using default_hash = murmur3_hash;

}
}
#endif // MURMUR3



#ifdef XXHASH
#include "hash/xx_hash.h"
#define HASHFCT dysect::hash::xx_hash
namespace dysect{
namespace hash {

    using default_hash = dysect::hash::xx_hash;

}
}
#endif // XXHASH
