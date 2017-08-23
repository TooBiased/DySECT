#pragma once

#include "include/config.h"

// cuckoo_simple tables

#ifdef CUCKOO_STANDARD
#define MULTI
#include "include/cuckoo_simple.h"
#define HASHTYPE cuckoo_standard
#endif // CUCKOO_STANDARD

#ifdef CUCKOO_STANDARD_INPLACE
#define MULTI
#include "include/cuckoo_simple.h"
#define HASHTYPE cuckoo_standard_inplace
#endif // CUCKOO_STANDARD_INPLACE



// cuckoo_deamortized table

#ifdef CUCKOO_DEAMORTIZED
#define MULTI
#include "include/cuckoo_deamortized.h"
#define HASHTYPE cuckoo_deamortized
#endif // CUCKOO_DEAMORTIZED



// cuckoo_dysect tables

#ifdef DYSECT
#define MULTI
#include "include/cuckoo_dysect.h"
#define HASHTYPE cuckoo_dysect
#endif // DYSECT

#ifdef DYSECT_INPLACE
#define MULTI
#include "include/cuckoo_dysect.h"
#define HASHTYPE cuckoo_dysect_inplace
#endif // DYSECT_INPLACE



// cuckoo_independent_2lvl table

#ifdef INDEPENDENT_2LVL
#define MULTI
#include "include/cuckoo_independent_2lvl.h"
#define HASHTYPE cuckoo_independent_2lvl
#endif // INDEPENDENT_2LVL



// cuckoo_overlap tables

#ifdef OVERLAP
#define MULTI
#include "include/cuckoo_overlap.h"
#define HASHTYPE cuckoo_overlap
#endif // OVERLAP

#ifdef OVERLAP_INPLACE
#define MULTI
#include "include/cuckoo_overlap.h"
#define HASHTYPE cuckoo_overlap_inplace
#endif // OVERLAP_INPLACE



// prob_hops tables

#ifdef HOPSCOTCH
#define HOPSCOTCH_CONFIG
#include "include/prob_hops.h"
#define  HASHTYPE prob_hopscotch
#endif // HOPSCOTCH

#ifdef HOPSCOTCH_INPLACE
#define HOPSCOTCH_CONFIG
#include "include/prob_hops.h"
#define  HASHTYPE prob_hopscotch_inplace
#endif // HOPSCOTCH_INPLACE



// prob_robin tables

#ifdef ROBIN
#define TRIV_CONFIG
#include "include/prob_robin.h"
#define  HASHTYPE prob_robin
#endif // ROBIN

#ifdef ROBIN_INPLACE
#define TRIV_CONFIG
#include "include/prob_robin.h"
#define  HASHTYPE prob_robin_inplace
#endif // ROBIN_INPLACE



// prob_simple tables

#ifdef LINEAR_DOUBLING
#define TRIV_CONFIG
#include "include/prob_simple.h"
#define  HASHTYPE prob_linear_doubling
#endif // LINEAR_DOUBLING

#ifdef LINEAR
#define TRIV_CONFIG
#include "include/prob_simple.h"
#define  HASHTYPE prob_linear
#endif // LINEAR

#ifdef LINEAR_INPLACE
#define TRIV_CONFIG
#include "include/prob_simple.h"
#define  HASHTYPE prob_linear_inplace
#endif // LINEAR_INPLACE



// NO TABLE CHOSEN!!! THEREFORE PRINT WARNING AND USE TEST TABLE!!!

#if (!MULTI       && \
     !TRIV_CONFIG && \
     !HOPSCOTCH_CONFIG)
#warning WARNING: No table chosen! Use
#define MULTI
#include "include/cuckoo_dysect.h"
#define  HASHTYPE cuckoo_dysect_inplace
#endif // NO TABLE IS DEFINED

#ifdef MULTI
#include "include/displacement_strategies/summary.h"
#endif // MULTI

#include "utils/commandline.h"

#ifdef NONCUCKOO
#define QUICK_MULTI
#endif  // NONCUCKOO


struct Chooser
{
#if defined QUICK_MULTI
    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    execute(CommandLine&, Types&& ... param)
    {
        Functor<Config<8,3,256,DisBFS> > f;
        return f(std::forward<Types>(param)...);
    }
#elif defined MULTI
    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    execute(CommandLine& c, Types&& ... param)
    {
        //return executeD<Functor, Hist, DisRWalk>     ( c, std::forward<Types>(param)...);
        ///*
        if      (c.boolArg("-bfs"))
            return executeD<Functor, Hist, DisBFS>     ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalk"))
            return executeD<Functor, Hist, DisRWalkOpt>   ( c, std::forward<Types>(param)...);
        else if (c.boolArg("-rwalkcyc"))
            return executeD<Functor, Hist, DisCycRWalk>( c, std::forward<Types>(param)...);

        std::cout << "ERROR: choose displacement Strategy (use triv)" << std::endl;
        return executeD<Functor, Hist, DisBFS>(c, std::forward<Types>(param)...);
        //*/
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
        case 1024:
            return executeDT<Functor, Hist, Displacer, 1024> (c, std::forward<Types>(param)...);
        // case 2048:
        //     return executeDT<Functor, Hist, Displacer, 2048> (c, std::forward<Types>(param)...);
        case 4096:
            return executeDT<Functor, Hist, Displacer, 4096> (c, std::forward<Types>(param)...);
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

#elif defined HOPSCOTCH_CONFIG
    template<template<class> class Functor, class Hist, class ... Types>
    inline static typename std::result_of<Functor<HopscotchConfig<> >(Types&& ...)>::type
    execute(CommandLine& c, Types&& ... param)
    {
        auto ns = c.intArg("-ns", HopscotchConfig<>::NeighborSize);
        switch (ns)
        {
        case 64:
            return executeN<Functor,64>(c, std::forward<Types>(param)...);
        case 128:
            return executeN<Functor,128>(c, std::forward<Types>(param)...);
        // case 8:
        //     return executeN<Functor, 8>(c, std::forward<Types>(param)...);
        // case 16:
        //     return executeN<Functor,16>(c, std::forward<Types>(param)...);
        // case 24:
        //     return executeN<Functor,24>(c, std::forward<Types>(param)...);
        // case 32:
        //     return executeN<Functor,32>(c, std::forward<Types>(param)...);
        // case 62:
        //     return executeN<Functor,62>(c, std::forward<Types>(param)...);
        default:
            constexpr auto tns = HopscotchConfig<>::NeighborSize;
            std::cout << "ERROR: unknown ns value (use "
                      << tns << ")" << std::endl;
            return executeN<Functor,tns>(c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, size_t NS, class ... Types>
    inline static typename std::result_of<Functor<HopscotchConfig<> >(Types&& ...)>::type
    executeN(CommandLine&
    #ifdef SPECIAL_HOPSCOTCH
             c
    #endif
             , Types&& ... param)
    {
        #ifdef SPECIAL_HOPSCOTCH
        double ratio = c.doubleArg("-rat", HopscotchConfig<>::GrowRatio_d);
        if (ratio < 1.101)
        {
            return executeNR<Functor, NS, std::ratio<11,10> >
                (std::forward<Types>(param)...);
        }
        else if (ratio < 1.1501)
        {
            return executeNR<Functor, NS, std::ratio<23,20> >
                (std::forward<Types>(param)...);
        }
        else if (ratio < 1.201)
        {
            return executeNR<Functor, NS, std::ratio<12,10> >
                (std::forward<Types>(param)...);
        }

        std::cout << "ERROR: unknown grow ratio (alpha) use "
                  << HopscotchConfig<>::GrowRatio_d << std::endl;
        return executeNR<Functor, NS, typename HopscotchConfig<>::GrowRatio>
            (std::forward<Types>(param)...);
        #else
        return executeNR<Functor, NS, std::ratio<1,2> >
            (std::forward<Types>(param)...);
        #endif
    }

    template<template<class> class Functor, size_t NS, class GRat, class ... Types>
    inline static typename std::result_of<Functor<HopscotchConfig<> >(Types&& ...)>::type
    executeNR(Types&& ... param)
    {
        Functor<HopscotchConfig<NS,GRat> > f;
        return f(std::forward<Types>(param)...);
    }

#elif defined TRIV_CONFIG
    template<template<class> class Functor, class Hist,
             class ... Types>
    inline static typename std::result_of<Functor<Config<> >(Types&& ...)>::type
    execute(CommandLine&, Types&& ... param)
    {
        Functor<TrivConfig<Hist> > f;
        return f(std::forward<Types>(param)...);
    }

#else

    template<template<class> class Functor, class Hist, class ... Types>
    inline static void execute(CommandLine&, Types&& ...)
    {
        std::cout << "some precompiler shit is broken" << std::endl;
    }
#endif
};
