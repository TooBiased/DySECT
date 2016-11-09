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
    using HashSplitter_t = typename Parent::HashSplitter_t;

    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    inline int insert(std::pair<Key, Data>, HashSplitter_t)
    {   return -1; }
};
