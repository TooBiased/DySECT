#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template <class Parent>
class dstrat_triv
{
public:
    using Key          = typename Parent::Key;
    using Data         = typename Parent::Data;
    using HashSplitter = typename Parent::HashSplitter;

    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    size_t  steps = 0;
    size_t* hist  = nullptr;

    bool insert(std::pair<Key, Data>, HashSplitter)
    {   return false;   }
};
