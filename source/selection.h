#pragma once

#include "include/config.h"

#if (!CUCKOO         && \
     !GROWS          && \
     !HOM2LVL        && \
     !TRIVGROW       && \
     !LINPROB        && \
     !SPACEPROB      && \
     !HOPSCOTCH      && \
     !SPACEHOPSCOTCH && \
     !MULTISIMPLE    && \
     !MULTIGROW      && \
     !MULTIHOM2LVL   && \
     !MULTITRIV      )
#define GROWS
#endif // NO HASHTYPE DEFINED => GROWS


#ifdef GROWS
#define SINGLE
#include "include/growing_cuckoo.h"
#define  HASHTYPE GrowingCuckoo
#endif // GROWS

#ifdef CUCKOO
#define SINGLE
#include "include/simple_cuckoo.h"
#define  HASHTYPE SimpleCuckoo
#endif // CUCKOO

#ifdef HOM2LVL
#define SINGLE
#include "include/hom_2lvl_cuckoo.h"
#define  HASHTYPE Hom2LvlCuckoo
#endif // HOM2LVL

#ifdef TRIVGROW
#define SINGLE
#include "include/triv_growing_cuckoo.h"
#define  HASHTYPE TrivGrowingCuckoo
#endif // TRIV_GROW

#ifdef LINPROB
#define NONCUCKOO
#include "include/lin_prob.h"
#define  HASHTYPE FastLinProb
#endif // LINPROB

#ifdef SPACEPROB
#define NONCUCKOO
#include "include/lin_prob.h"
#define  HASHTYPE SpaceLinProb
#endif

#ifdef HOPSCOTCH
#define NONCUCKOO
#include "include/hopscotch.h"
#define  HASHTYPE Hopscotch
#endif

#ifdef SPACEHOPSCOTCH
#define NONCUCKOO
#include "include/hopscotch.h"
#define  HASHTYPE SpaceHopscotch
#endif

#ifdef MULTISIMPLE
#define MULTI
#include "include/simple_multi_cuckoo.h"
#define HASHTYPE SimpleMultiCuckoo
#endif

#ifdef MULTIGROW
#define MULTI
#include "include/growing_multi_cuckoo.h"
#define HASHTYPE GrowingMultiCuckoo
#endif

#ifdef MULTIHOM2LVL
#define MULTI
#include "include/hom_2lvl_multi_cuckoo.h"
#define HASHTYPE Hom2LvlMultiCuckoo
#endif

#ifdef MULTITRIV
#define MULTI
#include "include/triv_growing_multi_cuckoo.h"
#define HASHTYPE TrivGrowingMultiCuckoo
#endif

#ifdef MULTI
#include "include/multistrategies/summary.h"
#endif // MULTI
#ifdef SINGLE
#include "include/strategies/summary.h"
#endif // SINGLE

#include "utils/commandline.h"

#ifndef NONCUCKOO
#define NORMAL
#endif  // NONCUCKOO

struct Chooser
{
#ifndef NORMAL
    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    execute(CommandLine&, Types&& ... param)
    {
        Functor<Config<> > f;
        return f(std::forward<Types>(param)...);
    }
#else

    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    execute(CommandLine& c, Types&& ... param)
    {
        if      (c.boolArg("-bfs"))
            return executeD<Functor, Hist, DisBFS>     ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalk"))
            return executeD<Functor, Hist, DisRWalk>   ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalkcyc"))
            return executeD<Functor, Hist, DisCycRWalk>( c, std::forward<Types>(param)...);

        std::cout << "ERROR: choose displacement Strategy (use triv)" << std::endl;
        return executeD<Functor, Hist, dstrat_triv>(c, std::forward<Types>(param)...);
    }

    template<template<class> class Functor, class Hist, template<class> class Displacer, class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    executeD(CommandLine& c, Types&& ... param)
    {
        auto tl = c.intArg("-tl", Config<>::tl);
        switch (tl)
        {
        // case 64:
        //     return executeDT<Functor, Hist, Displacer,   64> (c, std::forward<Types>(param)...);
        // case 128:
        //     return executeDT<Functor, Hist, Displacer,  128> (c, std::forward<Types>(param)...);
        // case 256:
        //     return executeDT<Functor, Hist, Displacer,  256> (c, std::forward<Types>(param)...);
        // case 512:
        //     return executeDT<Functor, Hist, Displacer,  512> (c, std::forward<Types>(param)...);
        // case 1024:
        //     return executeDT<Functor, Hist, Displacer, 1024> (c, std::forward<Types>(param)...);
        // case 2048:
        //     return executeDT<Functor, Hist, Displacer, 2048> (c, std::forward<Types>(param)...);
        default:
               constexpr auto ttl = Config<>::tl;
        //     std::cout << "ERROR: unknown TL value (use "
        //               << ttl << ")" << std::endl;
            return executeDT<Functor, Hist, Displacer, ttl>  (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist, template<class> class Displacer, size_t TL, class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    executeDT(CommandLine& c, Types&& ... param)
    {
        auto bs = c.intArg("-bs", Config<>::bs);
        switch (bs)
        {
        // case 4:
        //     return executeDTB<Functor, Hist, Displacer, TL,  4> (c, std::forward<Types>(param)...);
        // case 6:
        //     return executeDTB<Functor, Hist, Displacer, TL,  6> (c, std::forward<Types>(param)...);
        // case 8:
        //     return executeDTB<Functor, Hist, Displacer, TL,  8> (c, std::forward<Types>(param)...);
        // case 12:
        //     return executeDTB<Functor, Hist, Displacer, TL, 12> (c, std::forward<Types>(param)...);
        // case 16:
        //     return executeDTB<Functor, Hist, Displacer, TL, 16> (c, std::forward<Types>(param)...);
        default:
            constexpr auto tbs = Config<>::bs;
        //     std::cout << "ERROR: unknown BS value (use "
        //               << tbs << ")" << std::endl;
            return executeDTB<Functor, Hist, Displacer, TL, tbs> (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist,
             template<class> class Displacer, size_t TL, size_t BS,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    executeDTB(CommandLine&, Types&& ... param)
    {
        Functor<Config<BS,2,TL,Displacer,Hist> > f;
        return f(std::forward<Types>(param)...);
    }
#endif
};
