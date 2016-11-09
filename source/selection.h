#pragma once

#if (!CUCKOO && !GROWS && !HOM2LVL)
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

#ifdef HOM2LVL
#include "include/hom_2lvl_cuckoo.h"
#define HASHTYPE Hom2LvlCuckoo
#endif // HOM2LVL
