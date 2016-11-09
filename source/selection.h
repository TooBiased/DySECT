#pragma once

#if (!CUCKOO && !GROWS && !HOM2LVL)
#define GROWS
#endif // NO HASHTYPE DEFINED => GROWS


#ifdef GROWS
#include "include/growing_cuckoo.h"
#define HASHTYPE TGrowingCuckoo
#endif // GROWS

#ifdef CUCKOO
#include "include/simple_cuckoo.h"
#define HASHTYPE TSimpleCuckoo
#endif // CUCKOO

#ifdef HOM2LVL
#include "include/hom_2lvl_cuckoo.h"
#define HASHTYPE THom2LvlCuckoo
#endif // HOM2LVL
