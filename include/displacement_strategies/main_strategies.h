#pragma once

#include "dis_bfs1.h"
#include "dis_random_walk_optimistic.h"
#include "dis_trivial.h"

namespace dysect
{
namespace cuckoo_displacement
{

template <class c> using trivial     = dis_trivial<c>;
template <class c> using bfs         = dis_bfs1<c>;
template <class c> using random_walk = dis_random_walk_optimistic<c>;

} // namespace cuckoo_displacement
} // namespace dysect
