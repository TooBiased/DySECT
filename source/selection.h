#pragma once

#if (!CUCKOO && !GROWS && !HOM2LVL && !TRIV_GROW && !LINPROB)
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

#ifdef TRIV_GROW
#include "include/triv_growing_cuckoo.h"
#define HASHTYPE TrivGrowingCuckoo
#endif // TRIV_GROW

#ifdef LINPROB
#include "include/lin_prob.h"
#define HASHTYPE FastLinProbTable
#endif // LINPROB
