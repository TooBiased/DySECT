#pragma once

#include "dstrat_nbfs.h"
#include "dstrat_rwalk_anticycle.h"
#include "dstrat_rwalk_cyclic.h"
#include "dstrat_rwalk_optimistic.h"

template<class c> using DisBFS      = dstrat_nbfs<c>;
template<class c> using DisRWalk    = dstrat_rwalk_anticycle<c>;
template<class c> using DisCycRWalk = dstrat_rwalk_cyclic<c>;
template<class c> using DisRWalkOpt = dstrat_rwalk_optimistic<c>;
