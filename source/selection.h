#pragma once

#if (!CUCKOO && !GROWS)
#define GROWS
#endif // NO HASHTYPE DEFINED => GROWS


#ifdef GROWS
#include "include/growing_cuckoo.h"
#define HASHTYPE GrowingCuckoo
#endif // GROWS

#ifdef CUCKOO
#include "include/simple_cuckoo.h"
#define HASHTYPE SimpleCuckoo
#endif // CUCKOO

#ifdef CUCKOO_OLD
#include "include/cuckoo.h"
#define HASHTYPE CuckooTable
#endif // CUCKOO_OLD

#ifdef GROWS_OLD
#include "include/spacegrow.h"
#define HASHTYPE SpaceGrow
#endif // GROWS_OLD
