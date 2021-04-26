#pragma once

// cuckoo_simple tables

#ifdef MULTI_CUCKOO_STANDARD
#define MULTI
#include "include/cuckoo_simple.h"
#define HASHTYPE dysect::cuckoo_standard
#endif // CUCKOO_STANDARD

#ifdef MULTI_CUCKOO_STANDARD_INPLACE
#define MULTI
#include "include/cuckoo_simple.h"
#define HASHTYPE dysect::cuckoo_standard_inplace
#endif // CUCKOO_STANDARD_INPLACE



// cuckoo_deamortized table

#ifdef MULTI_CUCKOO_DEAMORTIZED
#define MULTI
#include "include/cuckoo_deamortized.h"
#define HASHTYPE dysect::cuckoo_deamortized
#endif // CUCKOO_DEAMORTIZED



// cuckoo_dysect tables

#ifdef MULTI_DYSECT
#define MULTI
#include "include/cuckoo_dysect.h"
#define HASHTYPE dysect::cuckoo_dysect
#endif // DYSECT

#ifdef MULTI_DYSECT_INPLACE
#define MULTI
#include "include/cuckoo_dysect.h"
#define HASHTYPE dysect::cuckoo_dysect_inplace
#endif // DYSECT_INPLACE



// cuckoo_independent_2lvl table

#ifdef MULTI_CUCKOO_INDEPENDENT_2LVL
#define MULTI
#include "include/cuckoo_independent_2lvl.h"
#define HASHTYPE dysect::cuckoo_independent_2lvl
#endif // INDEPENDENT_2LVL



// cuckoo_overlap tables

#ifdef MULTI_CUCKOO_OVERLAP
#define MULTI
#include "include/cuckoo_overlap.h"
#define HASHTYPE dysect::cuckoo_overlap
#endif // OVERLAP

#ifdef MULTI_CUCKOO_OVERLAP_INPLACE
#define MULTI
#include "include/cuckoo_overlap.h"
#define HASHTYPE dysect::cuckoo_overlap_inplace
#endif // OVERLAP_INPLACE



// prob_hops tables

#ifdef HOP_HOPSCOTCH
#define HOPSCOTCH_CONFIG
#include "include/prob_hops.h"
#define  HASHTYPE dysect::prob_hopscotch
#endif // HOPSCOTCH

#ifdef HOP_HOPSCOTCH_INPLACE
#define HOPSCOTCH_CONFIG
#include "include/prob_hops.h"
#define  HASHTYPE dysect::prob_hopscotch_inplace
#endif // HOPSCOTCH_INPLACE



// prob_robin tables

#ifdef TRIV_ROBIN
#define TRIV_CONFIG
#include "include/prob_robin.h"
#define  HASHTYPE dysect::prob_robin
#endif // ROBIN

#ifdef TRIV_ROBIN_INPLACE
#define TRIV_CONFIG
#include "include/prob_robin.h"
#define  HASHTYPE dysect::prob_robin_inplace
#endif // ROBIN_INPLACE



// prob_simple tables

// #ifdef LINEAR_DOUBLING
// #define TRIV_CONFIG
// #include "include/prob_simple.h"
// #define  HASHTYPE dysect::prob_linear_doubling
// #endif // LINEAR_DOUBLING

#ifdef TRIV_LINEAR
#define TRIV_CONFIG
#include "include/prob_simple.h"
#define  HASHTYPE dysect::prob_linear
#endif // LINEAR

#ifdef TRIV_LINEAR_INPLACE
#define TRIV_CONFIG
#include "include/prob_simple.h"
#define  HASHTYPE dysect::prob_linear_inplace
#endif // LINEAR_INPLACE

// prob_quadratic table

#ifdef TRIV_QUADRATIC
#define TRIV_CONFIG
#include "include/prob_quadratic.h"
#define  HASHTYPE dysect::prob_quadratic
#endif // QUADRATIC

#ifdef TRIV_QUADRATIC_INPLACE
#define TRIV_CONFIG
#include "include/prob_quadratic.h"
#define  HASHTYPE dysect::prob_quadratic_inplace
#endif // QUADRATIC

// multitable variants

#ifdef TRIV_MULTITABLE_LINEAR
#define TRIV_CONFIG
#include "include/prob_multitable_base.h"
#define HASHTYPE dysect::multitable_linear
#endif

#ifdef TRIV_MULTITABLE_ROBIN
#define TRIV_CONFIG
#include "include/prob_multitable_base.h"
#define HASHTYPE dysect::multitable_robin
#endif

#ifdef TRIV_MULTITABLE_QUADRATIC
#define TRIV_CONFIG
#include "include/prob_multitable_base.h"
#define HASHTYPE dysect::multitable_quadratic
#endif


// #ifdef MULTITABLE_QUADRATIC
// #define TRIV_CONFIG
// #include "include/prob_multitable_base.h"
// #define  HASHTYPE dysect::prob_quadratic
// #endif // QUADRATIC

// #ifdef MULTITABLE_CUCKOO
// #define MULTI
// #include "include/cuckoo_simple.h"
// #define HASHTYPE CuckooIndependentBase
// #endif


#ifdef TRIV_CHAINING
#define TRIV_CONFIG
#include "include/chaining.hpp"
#define HASHTYPE dysect::chaining
#endif




// NO TABLE CHOSEN!!! THEREFORE PRINT WARNING AND USE TEST TABLE!!!

/*
#if (!MULTI            &&                     \
     !TRIV_CONFIG      && \
     !HOPSCOTCH_CONFIG)

#warning WARNING: No table chosen! Use
#define MULTI
#include "include/cuckoo_dysect.h"
#define  HASHTYPE dysect::cuckoo_dysect_inplace
#endif // NO TABLE IS DEFINED
*/


#include "utils/command_line_parser.hpp"
namespace utm = utils_tm;

// #define QUICK_MULTI

struct Chooser
{
#if defined QUICK_MULTI
    template<template<class> class Functor, bool HistCount,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    execute(utm::command_line_parser&, Types&& ... param)
    {
        Functor<dysect::cuckoo_config<8,3,256,DisBFS> > f;
        return f(std::forward<Types>(param)...);
    }
#elif defined MULTI
    template<template<class> class Functor, bool HistCount,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    execute(utm::command_line_parser& c, Types&& ... param)
    {
        //return executeD<Functor, HistCount, DisRWalk>     ( c, std::forward<Types>(param)...);
        ///*
        if      (c.bool_arg("-bfs"))
            return executeD<Functor, HistCount, dysect::cuckoo_displacement::bfs>
                ( c, std::forward<Types>(param)...);
        else if (c.bool_arg("-rwalk"))
            return executeD<Functor, HistCount, dysect::cuckoo_displacement::random_walk>
                ( c, std::forward<Types>(param)...);

        std::cout << "ERROR: choose displacement Strategy (use triv)" << std::endl;
        return executeD<Functor, HistCount, dysect::cuckoo_displacement::trivial>
            (c, std::forward<Types>(param)...);
        //*/
    }

    template<template<class> class Functor, bool HistCount,
             template<class> class Displacer,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    executeD(utm::command_line_parser& c, Types&& ... param)
    {
        auto tl = c.int_arg("-tl", dysect::cuckoo_config<>::tl);
        switch (tl)
        {
        // case 64:
        //     return executeDT<Functor, HistCount, Displacer,   64> (c, std::forward<Types>(param)...);
        // case 128:
        //     return executeDT<Functor, HistCount, Displacer,  128> (c, std::forward<Types>(param)...);
        case 256:
            return executeDT<Functor, HistCount, Displacer,  256> (c, std::forward<Types>(param)...);
        // case 512:
        //     return executeDT<Functor, HistCount, Displacer,  512> (c, std::forward<Types>(param)...);
        // case 1024:
        //     return executeDT<Functor, HistCount, Displacer, 1024> (c, std::forward<Types>(param)...);
        // case 2048:
        //     return executeDT<Functor, HistCount, Displacer, 2048> (c, std::forward<Types>(param)...);
        // case 4096:
        //     return executeDT<Functor, HistCount, Displacer, 4096> (c, std::forward<Types>(param)...);
        default:
            constexpr auto ttl = dysect::cuckoo_config<>::tl;
            std::cout << "ERROR: unknown TL value (use "
                      << ttl << ")" << std::endl;
            return executeDT<Functor, HistCount, Displacer, ttl>  (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, bool HistCount,
             template<class> class Displacer, size_t TL,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    executeDT(utm::command_line_parser& c, Types&& ... param)
    {
        auto bs = c.int_arg("-bs", dysect::cuckoo_config<>::bs);
        switch (bs)
        {
        case 4:
            return executeDTB<Functor, HistCount, Displacer, TL,  4> (c, std::forward<Types>(param)...);
        // case 6:
        //     return executeDTB<Functor, HistCount, Displacer, TL,  6> (c, std::forward<Types>(param)...);
        case 8:
            return executeDTB<Functor, HistCount, Displacer, TL,  8> (c, std::forward<Types>(param)...);
        case 12:
            return executeDTB<Functor, HistCount, Displacer, TL, 12> (c, std::forward<Types>(param)...);
        // case 16:
        //     return executeDTB<Functor, HistCount, Displacer, TL, 16> (c, std::forward<Types>(param)...);
        default:
            constexpr auto tbs = dysect::cuckoo_config<>::bs;
            std::cout << "ERROR: unknown BS value (use "
                      << tbs << ")" << std::endl;
            return executeDTB<Functor, HistCount, Displacer, TL, tbs> (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, bool HistCount,
             template<class> class Displacer, size_t TL, size_t BS,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    executeDTB(utm::command_line_parser& c, Types&& ... param)
    {
        auto nh = c.int_arg("-nh", dysect::cuckoo_config<>::nh);
        switch (nh)
        {
        case 2:
            return executeDTBN<Functor, HistCount, Displacer, TL, BS, 2>
                (c, std::forward<Types>(param)...);
        case 3:
            return executeDTBN<Functor, HistCount, Displacer, TL, BS, 3>
                (c, std::forward<Types>(param)...);
        case 4:
            return executeDTBN<Functor, HistCount, Displacer, TL, BS, 4>
                (c, std::forward<Types>(param)...);
        default:
            constexpr auto tnh = dysect::cuckoo_config<>::nh;
            std::cout << "ERROR: unknown nh value (use "
                      << tnh << ")" << std::endl;
            return executeDTBN<Functor, HistCount, Displacer, TL, BS, tnh>
                (c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, bool HistCount,
             template<class> class Displacer, size_t TL, size_t BS, size_t NH,
             class ... Types>
    inline static typename std::result_of<Functor<dysect::cuckoo_config<> >(Types&& ...)>::type
    executeDTBN(utm::command_line_parser&, Types&& ... param)
    {
        if (HistCount)
        {
            Functor<dysect::cuckoo_config<BS,NH,TL,Displacer,
                                          false,
                                          dysect::hist_count> > f;
            return f(std::forward<Types>(param)...);
        }
        else
        {
            Functor<dysect::cuckoo_config<BS,NH,TL,Displacer,
                                          false,
                                          dysect::no_hist_count> > f;
            return f(std::forward<Types>(param)...);
        }
    }

#elif defined HOPSCOTCH_CONFIG
    template<template<class> class Functor, bool, class ... Types>
    inline static typename std::result_of<Functor<dysect::hopscotch_config<> >(Types&& ...)>::type
    execute(utm::command_line_parser& c, Types&& ... param)
    {
        auto ns = c.int_arg("-ns", dysect::hopscotch_config<>::neighborhood_size);
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
            constexpr auto tns = dysect::hopscotch_config<>::neighborhood_size;
            std::cout << "ERROR: unknown ns value (use "
                      << tns << ")" << std::endl;
            return executeN<Functor,tns>(c, std::forward<Types>(param)...);
        }
    }

    template<template<class> class Functor, size_t NS, class ... Types>
    inline static typename std::result_of<Functor<dysect::hopscotch_config<> >(Types&& ...)>::type
    executeN(utm::command_line_parser&, Types&& ... param)
    {
        Functor<dysect::hopscotch_config<NS> > f;
        return f(std::forward<Types>(param)...);
    }

#elif defined TRIV_CONFIG
    template<template<class> class Functor, bool, class ... Types>
    inline static typename std::result_of<Functor<dysect::triv_config>(Types&& ...)>::type
    execute(utm::command_line_parser&, Types&& ... param)
    {
        Functor<dysect::triv_config> f;
        return f(std::forward<Types>(param)...);
    }

#else

    template<template<class> class Functor, bool, class ... Types>
    inline static void execute(utm::command_line_parser&, Types&& ...)
    {
        std::cout << "some precompiler shit is broken" << std::endl;
    }
#endif
};
