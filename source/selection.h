#pragma once

#include "include/config.h"

#if (!CSIMPLE         && \
     !CEG2L           && \
     !CHOMOGENEOUS2L  && \
     !CINDEPENDENT2L  && \
     !LINPROB         && \
     !SPACEPROB       && \
     !HOPSCOTCH       && \
     !SPACEHOPSCOTCH)
#warning WARNING: No table chosen using! CEG2L
#define CEG2L
#endif // NO HASHTYPE DEFINED => GROWS

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

#ifdef CSIMPLE
#define MULTI
#include "include/cuckoo_simple.h"
#define HASHTYPE CuckooSimple
#endif

#ifdef CEG2L
#define MULTI
#include "include/cuckoo_eg2l.h"
#define HASHTYPE CuckooEG2L
#endif

#ifdef CHOMOGENEOUS2L
#define MULTI
#include "include/cuckoo_homogeneous2l.h"
#define HASHTYPE CuckooHomogeneous2L
#endif

#ifdef CINDEPENDENT2L
#define MULTI
#include "include/cuckoo_independent2l.h"
#define HASHTYPE CuckooIndependent2L
#endif

#ifdef MULTI
#include "include/displacement_strategies/summary.h"
#endif // MULTI

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
        case 256:
            return executeDT<Functor, Hist, Displacer,  256> (c, std::forward<Types>(param)...);
        // case 512:
        //     return executeDT<Functor, Hist, Displacer,  512> (c, std::forward<Types>(param)...);
        // case 1024:
        //     return executeDT<Functor, Hist, Displacer, 1024> (c, std::forward<Types>(param)...);
        // case 2048:
        //     return executeDT<Functor, Hist, Displacer, 2048> (c, std::forward<Types>(param)...);
        default:
            constexpr auto ttl = Config<>::tl;
            std::cout << "ERROR: unknown TL value (use "
                      << ttl << ")" << std::endl;
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
        case 4:
            return executeDTB<Functor, Hist, Displacer, TL,  4> (c, std::forward<Types>(param)...);
        // case 6:
        //     return executeDTB<Functor, Hist, Displacer, TL,  6> (c, std::forward<Types>(param)...);
        case 8:
            return executeDTB<Functor, Hist, Displacer, TL,  8> (c, std::forward<Types>(param)...);
        // case 12:
        //     return executeDTB<Functor, Hist, Displacer, TL, 12> (c, std::forward<Types>(param)...);
        // case 16:
        //     return executeDTB<Functor, Hist, Displacer, TL, 16> (c, std::forward<Types>(param)...);
        default:
            constexpr auto tbs = Config<>::bs;
            std::cout << "ERROR: unknown BS value (use "
                      << tbs << ")" << std::endl;
            return executeDTB<Functor, Hist, Displacer, TL, tbs> (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist,
             template<class> class Displacer, size_t TL, size_t BS,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    executeDTB(CommandLine& c, Types&& ... param)
    {
        auto nh = c.intArg("-nh", Config<>::nh);
        switch (nh)
        {
        case 2:
            return executeDTBN<Functor, Hist, Displacer, TL, BS, 2>
                (c, std::forward<Types>(param)...);
        case 3:
            return executeDTBN<Functor, Hist, Displacer, TL, BS, 3>
                (c, std::forward<Types>(param)...);
        // case 4:
        //     return executeDTBN<Functor, Hist, Displacer, TL, BS, 4>
        //         (c, std::forward<Types>(param)...);
        default:
            constexpr auto tnh = Config<>::nh;
            std::cout << "ERROR: unknown nh value (use "
                      << tnh << ")" << std::endl;
            return executeDTBN<Functor, Hist, Displacer, TL, BS, tnh>
                (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, class Hist,
             template<class> class Displacer, size_t TL, size_t BS, size_t NH,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    executeDTBN(CommandLine&, Types&& ... param)
    {
        Functor<Config<BS,NH,TL,Displacer,Hist> > f;
        return f(std::forward<Types>(param)...);
    }
#endif
};
