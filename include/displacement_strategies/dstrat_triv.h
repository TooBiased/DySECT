#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template <class Parent>
class dstrat_triv
{
private:
    using key_type      = typename Parent::key_type;
    using mapped_type   = typename Parent::mapped_type;
    using value_intern  = std::pair<key_type,mapped_type>;

    using Hashed_t      = typename Parent::Hashed_t;

public:
    dstrat_triv(Parent&, size_t, size_t) {}
    dstrat_triv(Parent&, dstrat_triv&&) {}

    inline std::pair<int, value_intern> insert(value_intern, Hashed_t)
    {   return std::make_pair(-1, nullptr); }
};
