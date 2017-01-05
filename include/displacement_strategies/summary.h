#pragma once

#include "dstrat_bfs.h"
#include "dstrat_rwalk_anticycle.h"
#include "dstrat_rwalk_cyclic.h"

template<class c> using DisBFS      = dstrat_bfs<c>;
template<class c> using DisRWalk    = dstrat_rwalk_anticycle<c>;
template<class c> using DisCycRWalk = dstrat_rwalk_cyclic<c>;
