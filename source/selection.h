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


#include "include/strategies/dstrat_bfs.h"
#include "include/strategies/dstrat_rwalk.h"
#include "include/strategies/dstrat_rwalk_cyclic.h"

#include "utils/commandline.h"

struct Chooser
{
    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<CuckooConfig<> >(Types&& ...)>::type
    execute(CommandLine& c, Types&& ... param)
    {
        if      (c.boolArg("-bfs"))
            return executeD<Functor, Hist, dstrat_bfs>         ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalk"))
            return executeD<Functor, Hist, dstrat_rwalk>       ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalkcyc"))
            return executeD<Functor, Hist, dstrat_rwalk_cyclic>( c, std::forward<Types>(param)...);

        std::cout << "ERROR: choose displacement Strategy (use triv)" << std::endl;
        return executeD<Functor, Hist, dstrat_triv>(c, std::forward<Types>(param)...);
    }

    template<template<class> class Functor, class Hist, template<class> class Displacer, class ... Types>
    inline static typename std::result_of<Functor<CuckooConfig<> >(Types&& ...)>::type
    executeD(CommandLine& c, Types&& ... param)
    {
        auto tl = c.intArg("-tl", CuckooConfig<>::tl);
        switch (tl)
        {
        case 64:
            return executeDT<Functor, Hist, Displacer, 64> (c, std::forward<Types>(param)...);
        case 256:
            return executeDT<Functor, Hist, Displacer, 256>  (c, std::forward<Types>(param)...);
        case 2048:
            return executeDT<Functor, Hist, Displacer, 2048> (c, std::forward<Types>(param)...);
        default:
            constexpr auto ttl = CuckooConfig<>::tl;
            std::cout << "ERROR: unknown TL value (use "
                      << ttl << ")" << std::endl;
            return executeDT<Functor, Hist, Displacer, ttl>  (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist, template<class> class Displacer, size_t TL, class ... Types>
    inline static typename std::result_of<Functor<CuckooConfig<> >(Types&& ...)>::type
    executeDT(CommandLine& c, Types&& ... param)
    {
        auto bs = c.intArg("-bs", CuckooConfig<>::bs);
        switch (bs)
        {
        case 4:
            return executeDTB<Functor, Hist, Displacer, TL, 4>  (c, std::forward<Types>(param)...);
        case 6:
            return executeDTB<Functor, Hist, Displacer, TL, 6>  (c, std::forward<Types>(param)...);
        case 8:
            return executeDTB<Functor, Hist, Displacer, TL, 8>  (c, std::forward<Types>(param)...);
        case 16:
            return executeDTB<Functor, Hist, Displacer, TL, 16> (c, std::forward<Types>(param)...);
        default:
            constexpr auto tbs = CuckooConfig<>::bs;
            std::cout << "ERROR: unknown BS value (use "
                      << tbs << ")" << std::endl;
            return executeDTB<Functor, Hist, Displacer, TL, tbs> (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist,
             template<class> class Displacer, size_t TL, size_t BS,
             class ... Types>
    inline static typename std::result_of<Functor<CuckooConfig<> >(Types&& ...)>::type
    executeDTB(CommandLine&, Types&& ... param)
    {
        Functor<CuckooConfig<BS,TL,Displacer,Hist> > f;
        return f(std::forward<Types>(param)...);
    }
};
