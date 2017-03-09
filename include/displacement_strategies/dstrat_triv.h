#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template <class Parent>
class dstrat_triv
{
public:
    using Key            = typename Parent::Key;
    using Data           = typename Parent::Data;
    using Pair_t         = std::pair<Key,Data>;
    using Hashed_t       = typename Parent::Hashed_t;

    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    inline std::pair<int, Pair_t> insert(Pair_t, Hashed_t)
    {   return std::make_pair(-1, nullptr); }
};
