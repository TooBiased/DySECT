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
    using Hashed_t       = typename Parent::Hashed_t;

    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    inline int insert(std::pair<Key, Data>, Hashed_t)
    {   return -1; }
};
