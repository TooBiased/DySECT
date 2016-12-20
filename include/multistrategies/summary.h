#pragma once

#include "dstrat_bfs.h"
#include "dstrat_rwalk.h"
#include "dstrat_rwalk_cyclic.h"

template<class c> using DisBFS      = dstrat_multi_bfs<c>;
template<class c> using DisRWalk    = dstrat_multi_rwalk<c>;
template<class c> using DisCycRWalk = dstrat_multi_rwalk_cyclic<c>;
